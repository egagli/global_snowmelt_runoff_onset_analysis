"""Per-tile UTM datacube construction: reproject a tile of the global
runoff-onset dataset to UTM 80 m and annotate it with the static ancillary
layers (DEM + slope/aspect, CHILI, seasonal snow class, ESA WorldCover,
forest cover fraction, GMBA range / HydroBASINS basin / continent IDs),
then tabulate (in memory) and emit the tile's PARTIAL SUMS for the
aggregation (gsro_analysis.aggregate.tile_partials) — the fleet product.

The ancillary layers are static per GRID generation — they do not change
with dataset version: built once per grid (compact integer encodings,
~8 MB/tile), joined with onset values per version. The per-tile pixel
table (parquet) is an OPT-IN product (process_tile(keep_pixels=True)):
the analyses never read it, the partials carry everything they need.
See pipeline/README.md.
"""

import numpy as np
import xarray as xr
import pystac_client
import xdem
import easysnowdata
import rasterio
import odc.stac
import planetary_computer
import logging
import geopandas as gpd
import rioxarray as rxr

from gsro_analysis import paths, settings


logger = logging.getLogger(__name__)


def setup_logging(level: int = logging.INFO) -> None:
    """
    Route root logging to the repo-root logs/ dir (never the caller's cwd —
    notebooks run from their own topic folders). OPT-IN since 2026-08-24:
    this used to run at module import, which meant merely importing
    gsro_analysis reconfigured the ROOT logger for the whole process — the
    same failure mode that once scattered gigabyte analysis.log files of
    azure/httpx chatter across the processing repo. Call it explicitly from
    entrypoints/notebooks that want file logging; library code just uses
    the module `logger`.
    """
    logging.basicConfig(filename=paths.logfile(), level=level,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def convert_water_year_dim_to_var(ds):
    for year in ds.water_year.values:
        ds[f'runoff_onset_WY{year}'] = ds['runoff_onset'].sel(water_year=year)

    ds = ds.drop_vars('runoff_onset').drop_vars('water_year')
    return ds

def add_topography(tile,ds):
    catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",modifier=planetary_computer.sign_inplace)
    search = catalog.search(collections=[f"cop-dem-glo-30"],intersects=tile.geobox.geographic_extent)
    dem_da = odc.stac.load(items=search.items(),like=ds,chunks={},resampling='bilinear')['data'].squeeze()
    dem_da = dem_da.rio.write_nodata(-32767,encoded=True).drop_vars('time') # compute for xdem stuff

    ds['dem'] = dem_da.compute()

    # [xDEM](https://xdem.readthedocs.io/en/stable/index.html) to calculate slope and aspect and topographic position index

    attributes = xdem.terrain.get_terrain_attribute(
        ds['dem'],
        resolution=ds['dem'].rio.resolution()[0],
        attribute=["aspect", "slope"], # , "topographic_position_index"
    )

    ds['aspect'] = xr.DataArray(attributes[0], dims=ds['dem'].dims, coords=ds['dem'].coords)
    ds['slope'] = xr.DataArray(attributes[1], dims=ds['dem'].dims, coords=ds['dem'].coords)
    # TPI? https://xdem.readthedocs.io/en/stable/gen_modules/xdem.DEM.topographic_position_index.html, https://tc.copernicus.org/articles/8/1989/2014/tc-8-1989-2014.pdf
    # maybe incorrect radius...
    #ds['tpi'] = xr.DataArray(attributes[2], dims=ds['dem'].dims, coords=ds['dem'].coords)

    # DAH?
    # alpha_max = 202.5 #only in northern hemisphere at specific latitude?
    # DAH_da = np.cos(np.deg2rad(alpha_max-aspect_da))*np.arctan(np.deg2rad(slope_da))

    return ds

CHILI_ASSET = "CSP/ERGo/1_0/Global/ALOS_CHILI"


def fetch_chili(bbox_gdf):
    """CHILI on its NATIVE ~90 m lattice cropped to the bbox, as a 0-1 index.

    xee >= 0.1 with an explicit pixel grid: ``easysnowdata.utils.get_ee_grid_params``
    snaps the bbox to the asset's own lattice, so the pixels are the asset's values,
    not a resample. Normalized by the asset's fixed 0-255 range — NOT the per-bbox
    min-max of ``easysnowdata.get_chili``, which made tiles incomparable (tile
    016_152's window spans 0-242, inflating its values by 5 %). Earth Engine must be
    initialized (settings.initialize_earthengine()).

    History (2026-09-03): easysnowdata <= 0.0.24 anchored the xee grid on the bbox
    corner, 0.46/0.51 px off the lattice, so every earlier CHILI layer was a
    nearest-neighbour resample of a half-pixel-shifted grid (max |diff| 0.42 vs the
    native values on 016_152). The v10 fleet builds from the native lattice."""
    import ee
    import xee  # noqa: F401  -- registers engine='ee'
    from easysnowdata.utils import get_ee_grid_params

    image = ee.Image(CHILI_ASSET)
    grid = get_ee_grid_params(image, bbox_gdf)
    da = xr.open_dataset(ee.ImageCollection(image), engine='ee', **grid)['constant'].isel(time=0, drop=True)
    da = (da.astype(np.float32) / np.float32(255.0))                        # dims (y, x), rioxarray-native names
    da = da.rio.write_crs(grid['crs']).rio.write_nodata(np.nan)
    da.attrs = {'long_name': 'Continuous Heat-Insolation Load Index, asset value / 255',
                'source': CHILI_ASSET, 'grid': 'native lattice',
                'data_citation': ('Theobald, D.M., Harrison-Atlas, D., Monahan, W.B., Albano, C.M. (2015). '
                                  'Ecologically-Relevant Maps of Landforms and Physiographic Diversity for '
                                  'Climate Adaptation Planning. PLoS ONE 10(12): e0143619.')}
    return da


def add_chili(tile, ds):
    chili_da = fetch_chili(tile.bbox_gdf)
    ds['chili'] = chili_da.rio.reproject_match(ds['dem'], resampling=rasterio.enums.Resampling.bilinear)
    return ds

def add_snow_class(tile,ds,mask_nodata=True):

    snow_classification = easysnowdata.remote_sensing.get_seasonal_snow_classification(tile.bbox_gdf,mask_nodata=True)
    ds['snow_classification'] = snow_classification.rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    return ds

def add_esa_worldcover(tile,ds,mask_nodata=True):

    esa_worldcover = easysnowdata.remote_sensing.get_esa_worldcover(tile.bbox_gdf, mask_nodata=True)
    ds['esa_worldcover'] = esa_worldcover.rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    return ds

def add_forest_cover(tile,ds,mask_nodata=True):

    #forest_cover_fraction = easysnowdata.remote_sensing.get_forest_cover_fraction(tile.bbox_gdf, mask_nodata=True)
    forest_cover_fraction = rxr.open_rasterio(
        settings.FOREST_COVER_FRACTION_URL,
        chunks=True,
        mask_and_scale=mask_nodata,
    ).squeeze().rio.clip_box(*tile.bbox_gdf.total_bounds, crs=tile.bbox_gdf.crs)
    ds['forest_cover_fraction'] = forest_cover_fraction.rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.bilinear)

    return ds

def add_mountain_range_and_basin_and_continent(tile,ds):
    from geocube.api.core import make_geocube

    # vector sources are downloaded ONCE into data/geometries/sources/ —
    # never fetched over the network per tile (BasinATLAS is ~1.9 GB and
    # figshare intermittently 202s, which corrupted a tile on 2026-08-24)
    gmba_zip = settings.cached_source(settings.GMBA_URL)
    gmba_clipped_gdf = gpd.read_file(
        gmba_zip, mask=tile.bbox_gdf).clip(tile.bbox_gdf)

    if gmba_clipped_gdf.empty:
        print(f"tile {tile.row},{tile.col} has no mountain ranges, filling with -9999")
        ds['GMBA_V2_ID'] = xr.full_like(ds['dem'], fill_value=-9999, dtype=np.int16)
    else:
        out_grid = make_geocube(
            vector_data=gmba_clipped_gdf,
            measurements=["GMBA_V2_ID"],
            resolution=(-0.0003, 0.0003),
        )

        ds['GMBA_V2_ID'] = out_grid['GMBA_V2_ID'].rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    basins_zip = settings.cached_source(settings.BASIN_ATLAS_URL,
                                        filename='BasinATLAS_Data_v10.gdb.zip',
                                        expected_md5=settings.BASIN_ATLAS_MD5)
    basins_gdf = gpd.read_file(basins_zip, mask=tile.bbox_gdf,
                               layer=settings.BASIN_ATLAS_LAYER)
    basins_clipped_gdf = basins_gdf.clip(tile.bbox_gdf)

    if basins_clipped_gdf.empty:
        print(f"tile {tile.row},{tile.col} has no basins, filling with -9999")
        ds['PFAF_ID'] = xr.full_like(ds['dem'], fill_value=-9999, dtype=np.int32)
    else:
        out_grid = make_geocube(
            vector_data=basins_clipped_gdf,
            measurements=["PFAF_ID"],
            resolution=(-0.0003, 0.0003),
        )
        
        ds['PFAF_ID'] = out_grid['PFAF_ID'].rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    continents_gdf = gpd.read_file(settings.cached_source(settings.CONTINENTS_URL))
    # NOTE: the integer continent encoding is the alphabetical order of
    # np.unique — aggregate.CONTINENTS_ENUM must match it, and the ancillary
    # build stamps it into the tile attrs so downstream code can assert.
    continents = list(np.unique(list(continents_gdf.CONTINENT)))
    categorical_enums = {'CONTINENT': continents}
    continents_clipped_gdf = continents_gdf.clip(tile.bbox_gdf)
    
    if continents_clipped_gdf.empty:
        print(f"tile {tile.row},{tile.col} has no continents, filling with -9999")
        ds['continent'] = xr.full_like(ds['dem'], fill_value=-9999, dtype=np.int16)
        
    else:
        out_grid = make_geocube(
            vector_data=continents_clipped_gdf,
            measurements=["CONTINENT"],
            resolution=(-0.0003, 0.0003),
            categorical_enums=categorical_enums
            
        ).where(lambda x: x >= 0)

        # geocube's "no category" code is -1 and it survives the mode resampling
        # as a value; it is nodata (48k coastal pixels on tile 016_152), not a class
        ds['continent'] = (out_grid['CONTINENT']
                           .rio.reproject_match(ds['dem'], resampling=rasterio.enums.Resampling.mode)
                           .where(lambda x: x >= 0))

    return ds


# ---------------------------------------------------------------------------
# Option A pipeline: static per-tile ancillary rasters + per-version tabulate
# (see pipeline/README.md). The ancillary is keyed by GRID generation and
# rebuilt only when the grid changes; tabulation joins onset values from any
# dataset version onto it by reproject_match, so alignment is exact by
# construction.

def tile_utm_template(config, row, col):
    """The tile's UTM 80 m target grid, derived WITHOUT onset data: zeros on
    the tile geobox (identical to a store clip) -> estimate_utm_crs ->
    rio.reproject(res=80) — the same derivation the v9-era datacube used."""
    import odc.geo.xr as odc_xr

    geobox = config.geobox_tiles[row, col]
    da = odc_xr.xr_zeros(geobox, dtype='float32')
    utm_crs = da.rio.estimate_utm_crs()
    return da.rio.reproject(utm_crs, resolution=80,
                            resampling=rasterio.enums.Resampling.bilinear)


def _utm_latlon_arrays(template):
    """Per-pixel WGS84 lat/lon of the UTM grid (exact transform of the pixel
    centers — supersedes the old bilinear-resampled original_lat/lon)."""
    import pyproj

    xx, yy = np.meshgrid(template.x.values, template.y.values)
    transformer = pyproj.Transformer.from_crs(template.rio.crs, "EPSG:4326",
                                              always_xy=True)
    lon, lat = transformer.transform(xx, yy)
    dims = ('y', 'x')
    return (xr.DataArray(lat.astype(np.float32), dims=dims),
            xr.DataArray(lon.astype(np.float32), dims=dims))


# Documented latitude coverage of each RASTER source, (south, north).
# OUTSIDE these bands a missing layer is expected and fills with nodata;
# INSIDE them a fetch failure is a real error and RAISES — never mask a
# transient failure as "no data" (a figshare 202 corrupted a mid-latitude
# tile on 2026-08-24 before this rule existed). Adjust after the
# far-north/-south dry-run tiles if a source proves narrower.
LAYER_LAT_COVERAGE = {
    'chili': (-70.0, 70.0),                  # per easysnowdata docstring
    'esa_worldcover': (-60.0, 82.75),
    'forest_cover_fraction': (-60.0, 78.25),  # PROBA-V LC100
}


def build_ancillary_tile(config, row, col, mask_nodata=True):
    """Static ancillary raster datacube for one tile on the UTM 80 m grid.

    Failure policy (fleet rule: failure = no output, never wrong output):
    a raster layer fills with nodata ONLY when the tile lies outside the
    source's documented latitude coverage (:data:`LAYER_LAT_COVERAGE`) or
    the source itself reports no data in bounds — any other failure
    (network, auth, API) raises so the tile is retried, not corrupted.
    The vector-ID layers (GMBA/PFAF/continent) handle legitimately-empty
    tiles internally and are never exception-wrapped. Attrs record the
    grid params, the continent enum, and provenance.
    """
    from rioxarray.exceptions import NoDataInBounds

    tile = config.get_tile(row, col)
    _, tile_south, _, tile_north = tile.bbox_gdf.total_bounds
    template = tile_utm_template(config, row, col)
    ds = template.to_dataset(name='_template')
    # lat/lon are NOT stored: they are an exact function of the UTM grid
    # (x, y, utm_crs) and are recomputed at tabulate time (_utm_latlon_arrays)

    def _fill(layers, reason):
        nonlocal ds
        logger.warning(f"tile {row},{col}: filling {layers} with nodata ({reason})")
        for name in layers:
            ds[name] = xr.full_like(template, np.nan)

    def _fetch(fn, *args, layers, coverage=None, **kwargs):
        nonlocal ds
        if coverage is not None and (tile_south > coverage[1]
                                     or tile_north < coverage[0]):
            _fill(layers, f"tile outside {fn.__name__} coverage {coverage}")
            return
        straddles = coverage is not None and (tile_south < coverage[0]
                                              or tile_north > coverage[1])
        try:
            ds = fn(*args, **kwargs)
        except NoDataInBounds as e:  # the source's own definitive empty signal
            _fill(layers, f"{fn.__name__}: {e}")
        except Exception as e:
            if straddles:  # partially past the coverage edge: expected
                _fill(layers, f"{fn.__name__} at coverage edge: {e}")
            else:          # inside documented coverage: a REAL failure
                raise

    _fetch(add_topography, tile, ds, layers=('dem', 'aspect', 'slope'))
    _fetch(add_chili, tile, ds, layers=('chili',),
           coverage=LAYER_LAT_COVERAGE['chili'])
    _fetch(add_snow_class, tile, ds, mask_nodata=mask_nodata,
           layers=('snow_classification',))
    _fetch(add_esa_worldcover, tile, ds, mask_nodata=mask_nodata,
           layers=('esa_worldcover',),
           coverage=LAYER_LAT_COVERAGE['esa_worldcover'])
    _fetch(add_forest_cover, tile, ds, mask_nodata=mask_nodata,
           layers=('forest_cover_fraction',),
           coverage=LAYER_LAT_COVERAGE['forest_cover_fraction'])
    # vector IDs: internal empty-handling only — exceptions must propagate
    ds = add_mountain_range_and_basin_and_continent(tile, ds)

    ds = ds.drop_vars('_template')
    # keep only serializable attrs: easysnowdata >= 0.0.25 attaches dict-valued attrs (snow-class
    # names/colours) that zarr would stringify and netCDF refuses
    for name in ds.data_vars:
        ds[name].attrs = {k: v for k, v in ds[name].attrs.items()
                          if isinstance(v, (str, bytes, int, float, list, tuple, np.ndarray, np.number))}
    continents_gdf = gpd.read_file(settings.cached_source(settings.CONTINENTS_URL))
    enum = list(np.unique(list(continents_gdf.CONTINENT)))
    ds.attrs.update({
        'tile_row': row, 'tile_col': col,
        'grid': f"{config.version}_grid",
        'utm_crs': str(template.rio.crs),
        'resolution_m': 80,
        'continent_enum': [f"{i}:{name}" for i, name in enumerate(enum)],
        'basin_atlas_layer': settings.BASIN_ATLAS_LAYER,   # the PFAF_ID level (6 since 2026-09-03)
    })
    try:
        from global_snowmelt_runoff_onset.provenance import collect_provenance
        ds.attrs['provenance'] = str(collect_provenance())
    except Exception:
        pass
    return ds


def ancillary_tile_path(config, row, col):
    """Azure path of a tile's ancillary zarr (container-relative)."""
    return (f"{settings.ANCILLARY_PREFIX}/{config.version}_grid/"
            f"tile_{row:03d}_{col:03d}.zarr")


def fresh_blob_fs(config):
    """See :func:`gsro_analysis.settings.fresh_blob_fs`."""
    return settings.fresh_blob_fs(config)


# Every variable a finished ancillary tile must contain. Layers outside a
# source's latitude coverage are still present — filled with nodata.
# (original_lat/lon were stored before 2026-08-26; tiles that still carry
# them are valid — the extra layers are ignored.)
ANCILLARY_VARS = ('dem', 'aspect', 'slope', 'chili', 'snow_classification',
                  'esa_worldcover', 'forest_cover_fraction', 'GMBA_V2_ID',
                  'PFAF_ID', 'continent')

# Compact on-disk encodings. The source layers arrive as float32/float64
# (make_geocube, easysnowdata, xdem) — 48 MB per tile as written before
# 2026-08-26, 8 MB with these (measured on tile 016_152; ~200 GB -> ~35 GB
# for the 4,320-tile grid). CF scale/offset + _FillValue mean the DECODED
# values are what tabulate_tile has always seen (chili is quantized to 1e-4,
# the precision the pixel table already rounds it to); nodata decodes to
# NaN and tabulate maps it to -9999 as before. One chunk per layer: tiles
# are always read whole, and 10 objects list/read faster than 220.
ANCILLARY_ENCODING = {
    'dem':                   {'dtype': 'int16', '_FillValue': -9999},
    'slope':                 {'dtype': 'uint8', '_FillValue': 255},
    'aspect':                {'dtype': 'int16', '_FillValue': -9999},
    'chili':                 {'dtype': 'uint16', '_FillValue': 65535,
                              'scale_factor': 1e-4, 'add_offset': 0.0},
    'snow_classification':   {'dtype': 'uint8', '_FillValue': 255},
    'esa_worldcover':        {'dtype': 'uint8', '_FillValue': 255},
    'forest_cover_fraction': {'dtype': 'uint8', '_FillValue': 255},
    'GMBA_V2_ID':            {'dtype': 'int32', '_FillValue': -9999},
    'PFAF_ID':               {'dtype': 'int32', '_FillValue': -9999},
    'continent':             {'dtype': 'int8', '_FillValue': -1},   # -1 = geocube's no-category code
}


def ancillary_encoding(ds):
    """zarr v3 encoding for save_ancillary_tile: the compact dtypes above,
    zstd, one chunk per layer. Integer-coded layers are rounded first (xarray
    rounds too, but be explicit; values already integral are unchanged)."""
    import zarr
    enc = {}
    comp = [zarr.codecs.BloscCodec(cname='zstd', clevel=5, shuffle='shuffle')]
    for name in ds.data_vars:
        e = dict(ANCILLARY_ENCODING.get(name, {}))
        if e and 'scale_factor' not in e:
            ds[name] = ds[name].round()
        # the source layers arrive with their own _FillValue/scale attrs and
        # encodings (easysnowdata, odc.stac); ours must win, so clear them
        for key in ('_FillValue', 'missing_value', 'scale_factor', 'add_offset', 'dtype'):
            ds[name].attrs.pop(key, None)
        ds[name].encoding = {}
        e['compressors'] = comp
        e['chunks'] = ds[name].shape
        enc[name] = e
    return enc

# Completion ledger: one marker blob per FINISHED tile, written only after
# to_zarr() returns. zarr writes its root metadata (zarr.json / .zmetadata)
# up front, so metadata presence can't distinguish a finished tile from one
# whose writer died mid-upload — tile 016_152 was exactly that on 2026-08-24
# (7 of 12 vars) and got skipped as "done". Existence checks must key on the
# ledger, never on zarr metadata. The ledger is a FLAT directory so the
# fleet dispatcher gets the whole done-set in one list call instead of
# probing 4,000+ tiles individually.

def ancillary_ledger_dir(config):
    return f"{settings.ANCILLARY_PREFIX}/{config.version}_grid/_complete"


def ancillary_marker_path(config, row, col):
    return f"{ancillary_ledger_dir(config)}/tile_{row:03d}_{col:03d}.json"


def ancillary_tile_complete(config, row, col):
    """True iff the tile's ancillary zarr finished writing (marker present)."""
    return config.azure_blob_fs.exists(ancillary_marker_path(config, row, col))


def completed_ancillary_tiles(config):
    """All (row, col) with a completion marker — one flat list call."""
    fs = config.azure_blob_fs
    ledger = ancillary_ledger_dir(config)
    if not fs.exists(ledger):
        return set()
    out = set()
    for p in fs.ls(ledger, detail=False):
        name = p.rsplit('/', 1)[-1]
        if name.startswith('tile_') and name.endswith('.json'):
            _, r, c = name[:-len('.json')].split('_')
            out.add((int(r), int(c)))
    return out


def _json_safe(value):
    """Coerce an attr value to something zarr's JSON attrs accept: numpy
    scalars/arrays -> Python, anything else unserializable -> str."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    return str(value)


def sanitize_attrs(ds):
    """The ancillary source layers (odc.stac, easysnowdata, geocube) attach
    attrs with numpy scalars / arbitrary objects that zarr's JSON attrs
    reject ('Invalid attribute in Dataset.attrs.') — coerce everything."""
    ds.attrs = {k: _json_safe(v) for k, v in ds.attrs.items()}
    for name in ds.variables:
        ds.variables[name].attrs = {
            k: _json_safe(v) for k, v in ds.variables[name].attrs.items()}
    return ds


def save_ancillary_tile(ds, config, row, col):
    missing = set(ANCILLARY_VARS) - set(ds.data_vars)
    if missing:
        raise ValueError(f"ancillary tile {row},{col} missing vars {sorted(missing)}")
    path = ancillary_tile_path(config, row, col)
    fs = config.azure_blob_fs
    marker = ancillary_marker_path(config, row, col)
    if fs.exists(marker):  # a rewrite: the old marker must not outlive the old zarr
        fs.rm(marker)
    if fs.exists(path):  # clear any partial previous attempt entirely
        fs.rm(path, recursive=True)
    fs.invalidate_cache()
    # consolidated=False on purpose: consolidation lists the group through
    # the (possibly stale) fsspec dircache — on 2026-08-25 that wrote a
    # 0-member consolidated_metadata and the store read back EMPTY. Thirteen
    # small arrays cost nothing to list at open time.
    ds = ds.drop_vars([v for v in ('original_lat', 'original_lon') if v in ds])
    ds = sanitize_attrs(ds)
    ds.to_zarr(fs.get_mapper(path), mode='w', consolidated=False,
               encoding=ancillary_encoding(ds))
    # verify-then-mark: the marker may only follow a successful re-read of
    # everything the tile must contain (fleet rule: no marker on failure).
    # In-process first (fresh fs), then a fresh interpreter — see
    # settings.verify_in_subprocess for why the in-process read can lie.
    if not (verify_and_mark_ancillary(config, row, col)
            or settings.verify_in_subprocess(config, 'datacube',
                                             'verify_and_mark_ancillary', row, col)):
        raise RuntimeError(f"tile {row},{col}: post-write verification failed "
                           f"({path})")
    logger.info(f"wrote {path}")
    return path


def verify_and_mark_ancillary(config, row, col):
    """Re-read the tile's zarr through a fresh fs; if every ancillary layer is
    present, write the completion marker and return True. Never raises."""
    import json
    path = ancillary_tile_path(config, row, col)
    try:
        fs = fresh_blob_fs(config)
        check = xr.open_zarr(fs.get_mapper(path), decode_coords='all',
                             chunks=None)
        missing = set(ANCILLARY_VARS) - set(check.data_vars)
        if missing:
            logger.warning(f"tile {row},{col}: store missing {sorted(missing)}")
            return False
        fs.pipe_file(ancillary_marker_path(config, row, col), json.dumps({
            'tile': [row, col], 'grid': f"{config.version}_grid",
            'data_vars': sorted(check.data_vars),
        }).encode())
        return True
    except Exception as e:  # noqa: BLE001
        logger.warning(f"tile {row},{col}: verification read failed: "
                       f"{type(e).__name__}: {e}")
        return False


def refresh_unit_layers(config, row, col, layers=('GMBA_V2_ID', 'PFAF_ID', 'continent')):
    """Stage 0b: re-rasterize the vector unit-id layers of a STORED ancillary
    tile and write them back in place (zarr ``r+``: the arrays are
    overwritten, nothing is deleted, the marker stays valid, the raster
    layers and their Earth Engine work are untouched). This is how a unit
    definition changes after the ancillary exists — e.g. the switch of
    ``PFAF_ID`` from HydroBASINS level 5 to level 6 (2026-09-03) on the
    dry-run tiles, or a future level on the whole grid (~30 s/tile, no EE).
    The tile's partials must be re-mapped afterwards (``process_tile``).
    Returns the refreshed layers as a Dataset."""
    tile = config.get_tile(row, col)
    stored = open_ancillary_tile(config, row, col)
    ds = stored[['dem']].copy()
    ds = add_mountain_range_and_basin_and_continent(tile, ds)
    new = ds[list(layers)]
    for name in layers:
        new[name] = new[name].round()
        for key in ('_FillValue', 'missing_value', 'scale_factor', 'add_offset', 'dtype'):
            new[name].attrs.pop(key, None)
        new[name].encoding = {}
        if name in ('GMBA_V2_ID', 'PFAF_ID'):   # nodata -> the stored fill, as tabulate expects
            new[name] = new[name].fillna(-9999)
    new = new.drop_vars([c for c in new.coords if c not in ('x', 'y')])
    fs = fresh_blob_fs(config)
    mapper = fs.get_mapper(ancillary_tile_path(config, row, col))
    new.to_zarr(mapper, mode='r+', consolidated=False)
    import zarr
    group = zarr.open_group(mapper, mode='r+')
    group.attrs['basin_atlas_layer'] = settings.BASIN_ATLAS_LAYER
    group.attrs['unit_layers_refreshed'] = ",".join(layers)
    fs.invalidate_cache()
    check = open_ancillary_tile(config, row, col)
    for name in layers:   # read back through the CF decoding tabulate uses
        got = check[name].values
        want = new[name].values
        if name in ('GMBA_V2_ID', 'PFAF_ID'):
            got = np.where(np.isfinite(got), got, -9999)
        if not np.array_equal(np.nan_to_num(got, nan=-1), np.nan_to_num(want, nan=-1)):
            raise RuntimeError(f"tile {row},{col}: {name} did not read back as written")
    logger.info(f"tile {row},{col}: refreshed {layers} in place")
    return check[list(layers)]


def open_ancillary_tile(config, row, col):
    # fresh fs on purpose: in a long-lived worker the shared fs can list
    # empty after the heavy build steps (see fresh_blob_fs)
    ds = xr.open_zarr(fresh_blob_fs(config).get_mapper(ancillary_tile_path(config, row, col)),
                      decode_coords='all', chunks=None)
    missing = set(ANCILLARY_VARS) - set(ds.data_vars)
    if missing:
        raise ValueError(
            f"ancillary tile {row},{col} is INCOMPLETE (missing {sorted(missing)}) "
            f"— a partial write; delete its zarr and marker and let the fleet rebuild it")
    ds = ds.drop_vars([v for v in ('original_lat', 'original_lon') if v in ds])
    return ds.rio.write_crs(ds.attrs['utm_crs'])


def tabulate_tile(config, row, col, ancillary_ds, global_ds=None):
    """Join one tile's runoff-onset values (any dataset version) onto its
    stored ancillary raster and tabulate.

    Reads with mask_and_scale=True and chunks=None (per-tile reads are
    small; CF decoding rescales runoff_onset_mad / temporal_resolution_median
    correctly and turns -9999 into NaN before the bilinear reprojection).
    Alignment is exact: reproject_match onto the ancillary grid, then an
    xr.merge(join='exact') tripwire.
    """
    from global_snowmelt_runoff_onset import store as gs_store

    if global_ds is None:
        global_ds = config.open_runoff_onset_dataset(chunks=None,
                                                     mask_and_scale=True)
    region = gs_store.tile_region_slices(config, row, col)
    tile_ds = global_ds.isel(region).drop_vars('temporal_resolution')
    tile_ds = tile_ds.astype(np.float32).compute()

    onset_utm = tile_ds.rio.reproject_match(
        ancillary_ds['dem'], resampling=rasterio.enums.Resampling.bilinear)
    onset_utm = convert_water_year_dim_to_var(onset_utm)
    ancillary_ds = ancillary_ds.copy()
    # exact WGS84 pixel-center coordinates of the stored UTM grid
    ancillary_ds['original_lat'], ancillary_ds['original_lon'] = \
        _utm_latlon_arrays(ancillary_ds['dem'])
    joined = xr.merge([onset_utm, ancillary_ds], join='exact',
                      combine_attrs='drop_conflicts')

    water_years = [int(y) for y in config.water_years]
    df = (joined.to_dataframe().reset_index()
          .dropna(subset='runoff_onset_median'))
    df = df.drop(columns=[c for c in ('spatial_ref',) if c in df])
    df['tile_row'] = row
    df['tile_col'] = col
    wy_cols = [f"runoff_onset_WY{y}" for y in water_years]
    df = df[["tile_row", "tile_col", "x", "y", "original_lat", "original_lon"]
            + wy_cols
            + ["runoff_onset_median", "temporal_resolution_median",
               "runoff_onset_mad", "dem", "aspect", "slope", "chili",
               "snow_classification", "esa_worldcover",
               "forest_cover_fraction", "GMBA_V2_ID", "PFAF_ID", "continent"]]

    int_cols = wy_cols + ["tile_row", "tile_col", "runoff_onset_median",
                          "dem", "aspect", "slope", "snow_classification",
                          "esa_worldcover", "forest_cover_fraction",
                          "GMBA_V2_ID", "continent"]
    df[int_cols] = (df[int_cols].replace([np.inf, -np.inf, np.nan], -9999)
                    .round().astype(np.int32))
    for c in ("tile_row", "tile_col", "runoff_onset_median", "dem", "aspect",
              "slope", "snow_classification", "esa_worldcover",
              "forest_cover_fraction", "continent"):
        df[c] = df[c].astype(np.int16)
    df["PFAF_ID"] = (df["PFAF_ID"].replace([np.inf, -np.inf, np.nan], -9999)
                     .astype(np.int32))
    for c in ("x", "y", "original_lat", "original_lon"):
        df[c] = df[c].astype(np.float32)
    df['chili'] = df['chili'].round(4).astype(np.float32)
    df['runoff_onset_mad'] = df['runoff_onset_mad'].round(2).astype(np.float32)
    df['temporal_resolution_median'] = (
        df['temporal_resolution_median'].round(2).astype(np.float32))
    for c in wy_cols:
        df[c] = df[c].astype(np.int16)
    return df


def parquet_tile_path(config, row, col):
    """Azure path of a tile's OPT-IN pixel table (container-relative). A
    parquet is a single blob (atomic on commit), so plain existence is a
    safe done-check."""
    return (f"{settings.ANALYSIS_PARQUET_PREFIX}/{config.version}/"
            f"tile_{row:03d}_{col:03d}.parquet")


def save_pixel_table(df, config, row, col):
    """Write a tile's pixel table: sorted by unit ids (row-group statistics
    then prune unit predicates), no pandas index column (it was 18 % of the
    file), zstd (-15 % vs snappy)."""
    path = parquet_tile_path(config, row, col)
    df = df.sort_values(['GMBA_V2_ID', 'PFAF_ID'], kind='stable')
    df.to_parquet(path, filesystem=config.azure_blob_fs, index=False,
                  compression='zstd', row_group_size=500_000)
    logger.info(f"wrote {path} ({len(df)} rows)")
    return path


def tabulate_and_save_tile(config, row, col, global_ds=None):
    """Tabulate one tile and write its pixel table (opt-in product)."""
    ancillary_ds = open_ancillary_tile(config, row, col)
    df = tabulate_tile(config, row, col, ancillary_ds, global_ds=global_ds)
    return save_pixel_table(df, config, row, col)


# ---------------------------------------------------------------------------
# the fleet product: per-tile partial sums (see aggregate.tile_partials)

def partials_tile_path(config, row, col):
    """Azure path of a tile's partial-sums parquet (single blob = atomic;
    existence = stage 1 done)."""
    return (f"{settings.PARTIALS_PREFIX}/{config.version}/"
            f"tile_{row:03d}_{col:03d}.parquet")


def completed_partials_tiles(config):
    """All (row, col) whose partials parquet exists — one flat list call."""
    fs = config.azure_blob_fs
    prefix = f"{settings.PARTIALS_PREFIX}/{config.version}"
    if not fs.exists(prefix):
        return set()
    out = set()
    for p in fs.ls(prefix, detail=False):
        name = p.rsplit('/', 1)[-1]
        if name.startswith('tile_') and name.endswith('.parquet'):
            _, r, c = name[:-len('.parquet')].split('_')
            out.add((int(r), int(c)))
    return out


def write_partials(df, config, row, col, filter_tags=None):
    """Partial sums of a tile's pixel table -> partials_tile_path() (single
    blob). Returns the number of partial rows."""
    from gsro_analysis import aggregate

    partials = aggregate.tile_partials(df, config.water_years, filter_tags=filter_tags)
    partials.insert(0, 'tile_col', np.int16(col))
    partials.insert(0, 'tile_row', np.int16(row))
    path = partials_tile_path(config, row, col)
    partials.to_parquet(path, filesystem=config.azure_blob_fs, index=False,
                        compression='zstd')
    logger.info(f"wrote {path} ({len(partials)} rows from {len(df)} pixels)")
    return len(partials)


def partials_from_pixel_table(config, row, col, filter_tags=None):
    """Backfill a tile's partials from an existing pixel parquet (no store
    read, no Earth Engine) — for tiles tabulated before the partials
    existed, e.g. the v10 dry-run set. The pixel table carries every column
    the map needs, so the result is identical to process_tile's."""
    import pandas as pd
    df = pd.read_parquet(parquet_tile_path(config, row, col),
                         filesystem=config.azure_blob_fs)
    return write_partials(df, config, row, col, filter_tags=filter_tags)


# ---------------------------------------------------------------------------
# starting a version over: the ONLY code in this repo that deletes Azure
# products (the fleet only ever adds). ERA5 stores and the dataset itself are
# never touched here.

VERSION_PRODUCTS = {
    'partials':      lambda c: f"{settings.PARTIALS_PREFIX}/{c.version}",           # fleet stage 1 (the aggregation input)
    'pixel_tables':  lambda c: f"{settings.ANALYSIS_PARQUET_PREFIX}/{c.version}",   # fleet, opt-in per-pixel parquets
    'ancillary_grid': lambda c: f"{settings.ANCILLARY_PREFIX}/{c.version}_grid",    # fleet stage 0 (per GRID; rebuilding = Earth Engine again)
    'aggregated_mirror': lambda c: f"{settings.AGGREGATED_PREFIX}/{c.version}",    # reduce --mirror copies of the cubes
    'era5_land':     lambda c: f"{settings.ERA5_LAND_PREFIX}/{c.version}",          # the ERA5-Land icechunk repo (acquisition + anomaly group;
                                                                                    # ERA5 Acquire rebuilds it, ~1 h) — NOT in the default reset
}
FLEET_PRODUCTS = VERSION_PRODUCTS  # historical alias


def version_products(config):
    """{name: (azure prefix, object count)} for every fleet product of this
    version/grid — a listing, nothing is touched."""
    fs = settings.fresh_blob_fs(config)
    out = {}
    for name, prefix_of in VERSION_PRODUCTS.items():
        prefix = prefix_of(config)
        out[name] = (prefix, len(fs.find(prefix)) if fs.exists(prefix) else 0)
    return out


def reset_version(config, what=('partials', 'pixel_tables', 'ancillary_grid'),
                  confirm=False):
    """Delete the chosen fleet products of this version on Azure so the next
    dispatch rebuilds every tile from scratch (the dispatcher lists remaining
    work from these prefixes). Prints the plan and does NOTHING unless
    ``confirm=True``. ``'era5_land'`` must be asked for explicitly (then the
    ERA5 Acquire workflow re-acquires every water year; its start_fresh box does the same). Never deletes the icechunk dataset,
    the pyramid, or anything outside VERSION_PRODUCTS."""
    unknown = set(what) - set(VERSION_PRODUCTS)
    if unknown:
        raise ValueError(f"unknown products {sorted(unknown)}; choose from {list(VERSION_PRODUCTS)}")
    fs = settings.fresh_blob_fs(config)
    plan = {k: v for k, v in version_products(config).items() if k in what}
    for name, (prefix, n) in plan.items():
        print(f"{'DELETE' if confirm else 'would delete'} {name:18s} {prefix}  ({n} objects)")
    if not confirm:
        print("dry run — pass confirm=True to delete")
        return plan
    for name, (prefix, n) in plan.items():
        if n:
            fs.rm(prefix, recursive=True)
    fs.invalidate_cache()
    config.azure_blob_fs.invalidate_cache()
    print("done; the next fleet dispatch rebuilds every tile in the work list")
    return plan


def process_tile(config, row, col, global_ds=None, keep_pixels=False,
                 filter_tags=None):
    """The whole per-tile job: ancillary (built only if its completion
    marker is missing) -> in-memory pixel table -> partial sums written to
    partials_tile_path(); the pixel table itself is written only with
    ``keep_pixels``. Failure = exception = no partials blob (the fleet rule).
    Returns (n_pixels, n_partial_rows)."""
    if not ancillary_tile_complete(config, row, col):
        ds = build_ancillary_tile(config, row, col)
        save_ancillary_tile(ds, config, row, col)
    ancillary_ds = open_ancillary_tile(config, row, col)
    df = tabulate_tile(config, row, col, ancillary_ds, global_ds=global_ds)
    if keep_pixels:
        save_pixel_table(df, config, row, col)
    n_rows = write_partials(df, config, row, col, filter_tags=filter_tags)
    return len(df), n_rows


def open_coarse_onset(config, level=5):
    """Coarse-resolution runoff-onset dataset for masking/joining against
    other gridded products: the public multiscale pyramid (level n ~
    80 m * 2**n; 5 ~ 2.6 km, 7 ~ 10 km ~ ERA5-Land), anonymous."""
    from global_snowmelt_runoff_onset.pyramid import open_pyramid_level
    return open_pyramid_level(config, level).rio.write_crs("EPSG:4326")
