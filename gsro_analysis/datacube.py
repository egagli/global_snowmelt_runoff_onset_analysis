"""Per-tile ancillary and tabulation engine of the analysis pipeline, in two waves
(pipeline/README.md):

**Get ancillary data** (once per GRID generation): the eight source layers
(:data:`ANCILLARY_LAYERS`) resampled onto the tile's window of the DATASET grid
(EPSG:4326, 0.00072 deg, 2048 x 2048 per tile) and written into ONE icechunk
repository on the global dataset geobox with compact integer encodings
(:data:`ANCILLARY_ENCODING`; one chunk per tile per layer). One commit per tile
with metadata is the ledger (gsro_analysis.ledger, the production fleet
pattern): :func:`initialize_ancillary_store`, :func:`build_ancillary_window`,
:func:`write_ancillary_tile`, :func:`completed_ancillary_tiles`,
:func:`open_ancillary_window`.

**Process tiles to parquets** (per dataset version): :func:`process_tile` reads
the tile's window from the dataset store and from the ancillary store (the same
pixel indices: both live on the dataset grid), reprojects BOTH onto the tile's
80 m UTM grid (:func:`tile_utm_template`; onset, DEM, CHILI and forest cover
bilinear, the categorical and id layers nearest), derives slope and aspect from
the UTM DEM (xdem), tabulates the pixels (:func:`tabulate_tile`) and writes the
tile's partial sums (:func:`write_partials`; the pixel table is opt-in). Every
row is a ~6,400 m2 UTM pixel, so pixel counts stay area and no area weight is
needed; the price is that the continuous layers are resampled twice.

Unit-definition changes after the ancillary exists are :func:`refresh_unit_layers`
(re-rasterize the id layers into the store, one commit, no Earth Engine) plus a
re-map. :func:`reset_version` is the only deleting code (dry run by default).
"""

import logging
import time

import easysnowdata
import geopandas as gpd
import numpy as np
import odc.stac
import planetary_computer
import pystac_client
import rasterio
import rioxarray as rxr
import xarray as xr
import xdem

from gsro_analysis import ledger, paths, settings

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


# ---------------------------------------------------------------------------
# the source layers, each resampled onto the target grid carried by ds['dem'] (or the template)

def add_dem(tile, ds):
    """Copernicus DEM GLO-30 (Planetary Computer) resampled bilinearly onto the target grid
    (``like=ds``). Slope and aspect are NOT derived here: they need a metric grid and are
    computed at process time from the UTM-reprojected DEM (:func:`terrain_derivatives`)."""
    catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1",
                                        modifier=planetary_computer.sign_inplace)
    search = catalog.search(collections=["cop-dem-glo-30"], intersects=tile.geobox.geographic_extent)
    items = list(search.items())
    if not items:
        # the source's own definitive no-data signal (build_ancillary_window fills the layer with nodata):
        # the public GLO-30 release has no tiles over Armenia and Azerbaijan (work-list tiles 29,152 /
        # 29,153 / 29,154 on 2026-09-04) — such a tile contributes no elevation-binned statistics
        from rioxarray.exceptions import NoDataInBounds
        raise NoDataInBounds("Copernicus DEM GLO-30: no STAC items intersect the tile "
                             "(open water, or the Armenia/Azerbaijan exclusion of the public release)")
    dem_da = odc.stac.load(items=items, like=ds, chunks={}, resampling='bilinear')['data'].squeeze()
    dem_da = dem_da.rio.write_nodata(-32767, encoded=True).drop_vars('time')
    ds['dem'] = dem_da.compute()
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
    # a static asset carries no system:time_start; stamp one so xee builds its (dropped) time axis quietly
    collection = ee.ImageCollection([image.set('system:time_start', 0)])
    da = xr.open_dataset(collection, engine='ee', **grid)['constant'].isel(time=0, drop=True)
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


def add_mountain_range_and_basin_and_continent(tile, ds, log=print):
    from geocube.api.core import make_geocube

    # vector sources are downloaded ONCE into analyses/<unit>/data/geometries/ —
    # never fetched over the network per tile (BasinATLAS is ~2.7 GB and
    # figshare intermittently 202s, which corrupted a tile on 2026-08-24)
    gmba_clipped_gdf = gpd.read_file(
        settings.gmba_zip(), mask=tile.bbox_gdf).clip(tile.bbox_gdf)

    if gmba_clipped_gdf.empty:
        log(f"tile {tile.row},{tile.col}: no mountain ranges in the tile, id layer filled with -9999")
        ds['GMBA_V2_ID'] = xr.full_like(ds['dem'], fill_value=-9999, dtype=np.int16)
    else:
        out_grid = make_geocube(
            vector_data=gmba_clipped_gdf,
            measurements=["GMBA_V2_ID"],
            resolution=(-0.0003, 0.0003),
        )

        ds['GMBA_V2_ID'] = out_grid['GMBA_V2_ID'].rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    basins_gdf = gpd.read_file(settings.basin_atlas_gdb(), mask=tile.bbox_gdf,
                               layer=settings.BASIN_ATLAS_LAYER)
    basins_clipped_gdf = basins_gdf.clip(tile.bbox_gdf)

    if basins_clipped_gdf.empty:
        log(f"tile {tile.row},{tile.col}: no basins in the tile, id layer filled with -9999")
        ds['PFAF_ID'] = xr.full_like(ds['dem'], fill_value=-9999, dtype=np.int32)
    else:
        out_grid = make_geocube(
            vector_data=basins_clipped_gdf,
            measurements=["PFAF_ID"],
            resolution=(-0.0003, 0.0003),
        )
        
        ds['PFAF_ID'] = out_grid['PFAF_ID'].rio.reproject_match(ds['dem'],resampling=rasterio.enums.Resampling.mode)

    continents_gdf = gpd.read_file(settings.continents_zip())
    # NOTE: the integer continent encoding is the alphabetical order of
    # np.unique — aggregate.CONTINENTS_ENUM must match it, and the ancillary
    # build stamps it into the tile attrs so downstream code can assert.
    continents = list(np.unique(list(continents_gdf.CONTINENT)))
    categorical_enums = {'CONTINENT': continents}
    continents_clipped_gdf = continents_gdf.clip(tile.bbox_gdf)
    
    if continents_clipped_gdf.empty:
        log(f"tile {tile.row},{tile.col}: no continents in the tile, id layer filled with -9999")
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
# the tile's UTM 80 m grid: the TABULATION grid (the rows of the pixel table are its pixels)

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


# ---------------------------------------------------------------------------
# the ancillary store: one icechunk repository per grid generation on the dataset grid

# the eight layers fetched from their sources and stored on the dataset grid
ANCILLARY_LAYERS = ('dem', 'chili', 'snow_classification', 'esa_worldcover',
                    'forest_cover_fraction', 'GMBA_V2_ID', 'PFAF_ID', 'continent')
# derived at process time from the UTM-reprojected DEM (they need a metric grid)
DERIVED_LAYERS = ('aspect', 'slope')
# every ancillary column of the pixel table
ANCILLARY_VARS = ('dem', 'aspect', 'slope', 'chili', 'snow_classification',
                  'esa_worldcover', 'forest_cover_fraction', 'GMBA_V2_ID',
                  'PFAF_ID', 'continent')

# Compact on-disk encodings (CF scale/offset + _FillValue, so the DECODED values are what
# tabulate_tile sees: chili quantized to 1e-4, nodata -> NaN). ~56 MB raw per tile.
ANCILLARY_ENCODING = {
    'dem':                   {'dtype': 'int16', '_FillValue': -9999},
    'chili':                 {'dtype': 'uint16', '_FillValue': 65535,
                              'scale_factor': 1e-4, 'add_offset': 0.0},
    'snow_classification':   {'dtype': 'uint8', '_FillValue': 255},
    'esa_worldcover':        {'dtype': 'uint8', '_FillValue': 255},
    'forest_cover_fraction': {'dtype': 'uint8', '_FillValue': 255},
    'GMBA_V2_ID':            {'dtype': 'int32', '_FillValue': -9999},
    'PFAF_ID':               {'dtype': 'int32', '_FillValue': -9999},
    'continent':             {'dtype': 'int8', '_FillValue': -1},   # -1 = geocube's no-category code
}
LAYER_ATTRS = {
    'dem': {'long_name': 'Copernicus DEM GLO-30 elevation', 'units': 'm'},
    'chili': {'long_name': 'Continuous Heat-Insolation Load Index, asset value / 255',
              'source': 'CSP/ERGo/1_0/Global/ALOS_CHILI'},
    'snow_classification': {'long_name': 'Sturm & Liston 2021 seasonal snow classification'},
    'esa_worldcover': {'long_name': 'ESA WorldCover v200 land cover class'},
    'forest_cover_fraction': {'long_name': 'PROBA-V LC100 tree cover fraction', 'units': 'percent'},
    'GMBA_V2_ID': {'long_name': 'GMBA mountain inventory v2 (standard 300) range id'},
    'PFAF_ID': {'long_name': 'HydroBASINS Pfafstetter id (BasinATLAS)'},
    'continent': {'long_name': 'continent index (see continent_enum)'},
}
# how each stored layer moves from the dataset grid onto the 80 m UTM grid at process time
UTM_RESAMPLING = {
    'dem': rasterio.enums.Resampling.bilinear,
    'chili': rasterio.enums.Resampling.bilinear,
    'forest_cover_fraction': rasterio.enums.Resampling.bilinear,
    'snow_classification': rasterio.enums.Resampling.nearest,
    'esa_worldcover': rasterio.enums.Resampling.nearest,
    'GMBA_V2_ID': rasterio.enums.Resampling.nearest,
    'PFAF_ID': rasterio.enums.Resampling.nearest,
    'continent': rasterio.enums.Resampling.nearest,
}
KIND_ANCILLARY_TILE = 'ancillary_tile'


def ancillary_repo_prefix(config):
    """Container-qualified Azure prefix of the grid generation's ancillary icechunk repository."""
    return f"{settings.ANCILLARY_PREFIX}/{config.version}_grid/ancillary"


def ancillary_repo_exists(config, local_store=None):
    return ledger.repo_exists(config, ancillary_repo_prefix(config), local_store)


def open_ancillary_repo(config, local_store=None):
    return ledger.open_repo(config, ancillary_repo_prefix(config), local_store, what='ancillary repository')


def tile_geo_template(config, row, col):
    """The tile's window of the dataset grid (EPSG:4326, 0.00072 deg, 2048 x 2048 except edge
    tiles) as a float32 zeros DataArray with CRS: the build target of the ancillary layers."""
    import odc.geo.xr as odc_xr
    da = odc_xr.xr_zeros(config.geobox_tiles[row, col], dtype='float32')
    if 'y' in da.dims:
        da = da.rename({'y': 'latitude', 'x': 'longitude'})
    return da.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude').rio.write_crs('EPSG:4326')


def build_ancillary_template(config):
    """Lazy all-NaN template of the eight layers on the GLOBAL dataset geobox plus its Zarr v3
    encoding (compact ints, CF fill/scale, zstd, one chunk per tile per layer; the zarr
    fill_value is the encoded nodata so never-written tiles decode to NaN). Written
    metadata-only by :func:`initialize_ancillary_store`."""
    import dask.array as da
    import odc.geo.xr as odc_xr
    import zarr
    geobox = config.global_geobox
    tile_dim = config.spatial_chunk_dim_zarr_output
    layers = {}
    for name in ANCILLARY_LAYERS:
        arr = odc_xr.wrap_xr(da.full(geobox.shape.yx, np.nan, dtype='float32', chunks=(tile_dim, tile_dim)), geobox)
        layers[name] = arr.rename(name)
    ds = xr.Dataset(layers)
    if 'y' in ds.dims:
        ds = ds.rename({'y': 'latitude', 'x': 'longitude'})
    for name in ANCILLARY_LAYERS:
        ds[name].attrs.update(LAYER_ATTRS[name])
        ds[name].attrs['grid_mapping'] = 'spatial_ref'
    continents_gdf = gpd.read_file(settings.continents_zip())
    enum = list(np.unique(list(continents_gdf.CONTINENT)))
    ds.attrs.update({
        'title': 'Static ancillary layers of the global snowmelt runoff onset analyses, on the dataset grid',
        'grid': f"{config.version}_grid", 'crs': 'EPSG:4326', 'resolution_deg': float(config.resolution),
        'tile_dim': int(tile_dim), 'layers': list(ANCILLARY_LAYERS),
        'derived_at_process_time': list(DERIVED_LAYERS),
        'continent_enum': [f"{i}:{name}" for i, name in enumerate(enum)],
        'basin_atlas_layer': settings.BASIN_ATLAS_LAYER,
        'chili_normalization': 'asset value / 255 (fixed 0-255 range, native lattice)',
    })
    comp = [zarr.codecs.BloscCodec(cname='zstd', clevel=5, shuffle='shuffle')]
    encoding = {}
    for name in ANCILLARY_LAYERS:
        e = dict(ANCILLARY_ENCODING[name])
        encoding[name] = {**e, 'fill_value': e['_FillValue'], 'chunks': (tile_dim, tile_dim), 'compressors': comp}
    return sanitize_attrs(ds), encoding


def initialize_ancillary_store(config, start_fresh=False, local_store=None, log=print):
    """Return the grid generation's ancillary repository, creating it with the empty template
    if it does not exist. ``start_fresh=True`` DELETES an existing repository first (the
    'Get ancillary data' workflow's off-by-default box)."""
    prefix = ancillary_repo_prefix(config)
    if ledger.repo_exists(config, prefix, local_store):
        if not start_fresh:
            log(f"ancillary repository exists: {local_store or prefix}")
            return ledger.open_repo(config, prefix, local_store)
        ledger.delete_repo(config, prefix, local_store, log=log)
    repo = ledger.create_repo(config, prefix, local_store)
    template, encoding = build_ancillary_template(config)
    session = repo.writable_session(ledger.BRANCH)
    template.to_zarr(session.store, mode='w', zarr_format=3, compute=False, write_empty_chunks=False,
                     consolidated=False, encoding=encoding)
    session.commit(f"initialize empty ancillary store, {config.version}_grid",
                   metadata={'schema': ledger.SCHEMA, 'kind': 'init', 'grid': f"{config.version}_grid",
                             'layers': list(ANCILLARY_LAYERS), 'tile_dim': int(config.spatial_chunk_dim_zarr_output),
                             'shape': [int(v) for v in config.global_geobox.shape.yx],
                             'provenance': ledger.provenance()})
    log(f"initialized {local_store or prefix}: empty template {tuple(config.global_geobox.shape.yx)}, "
        f"{len(ANCILLARY_LAYERS)} layers, one chunk per tile")
    return repo


def ancillary_records(repo):
    """Newest -> oldest ancillary-tile commit records (kind ``ancillary_tile``)."""
    return [r for r in ledger.commit_records(repo) if r.get('kind') == KIND_ANCILLARY_TILE]


def completed_ancillary_tiles(config, repo=None, local_store=None):
    """{(row, col)} with an ancillary commit: the fold over the store's commit history
    (newest wins; a refreshed tile stays complete). Empty when the repository does not exist."""
    if repo is None:
        if not ancillary_repo_exists(config, local_store):
            return set()
        repo = open_ancillary_repo(config, local_store)
    return {tuple(int(v) for v in r['tile']) for r in ancillary_records(repo)}


def ancillary_tile_complete(config, row, col, repo=None):
    return (row, col) in completed_ancillary_tiles(config, repo)


def _with_transient_retries(fn, *args, log=None, tries=3, delays=(20, 60), **kwargs):
    """Call a layer fetch, retrying only on the transient failures of the remote catalogs
    (a Planetary Computer STAC gateway error page, a dropped connection): two of six wave-2
    batches lost one tile each to a ~03:00 UTC blip on 2026-09-04. Anything else raises at
    once; after ``tries`` the transient error raises too (failure = no commit)."""
    import requests
    import pystac_client.exceptions
    transient = (pystac_client.exceptions.APIError, requests.exceptions.RequestException,
                 ConnectionError, TimeoutError)
    for attempt in range(tries):
        try:
            return fn(*args, **kwargs)
        except transient as e:
            if attempt == tries - 1:
                raise
            delay = delays[min(attempt, len(delays) - 1)]
            (log or logger.warning)(f"{fn.__name__}: transient {type(e).__name__}, retry {attempt + 2}/{tries} in {delay}s")
            time.sleep(delay)


def build_ancillary_window(config, row, col, mask_nodata=True, log=None):
    """The tile's eight source layers on its window of the DATASET grid (wave 'Get ancillary
    data'). Failure policy (fleet rule: failure = no output, never wrong output): a raster
    layer fills with nodata ONLY when the tile lies outside the source's documented latitude
    coverage (:data:`LAYER_LAT_COVERAGE`) or the source itself reports no data in bounds;
    any other failure (network, auth, API) raises so the tile is retried, not corrupted.
    The vector-id layers handle legitimately-empty tiles internally and are never
    exception-wrapped. ``log`` receives one line per step."""
    from rioxarray.exceptions import NoDataInBounds
    log = log or logger.info
    tile = config.get_tile(row, col)
    _, tile_south, _, tile_north = tile.bbox_gdf.total_bounds
    template = tile_geo_template(config, row, col)
    ds = template.to_dataset(name='_template')
    log(f"tile {row},{col}: ancillary target = the tile's window of the dataset grid "
        f"(EPSG:4326, {template.sizes['latitude']} x {template.sizes['longitude']} px at 0.00072 deg)")
    step_names = {
        'add_dem': 'Copernicus DEM GLO-30 -> bilinear',
        'add_chili': 'CHILI (Earth Engine, native lattice via xee) -> bilinear',
        'add_snow_class': 'seasonal snow classification -> mode',
        'add_esa_worldcover': 'ESA WorldCover -> mode',
        'add_forest_cover': 'forest cover fraction (PROBA-V LC100) -> bilinear',
    }

    def _fill(layers, reason):
        nonlocal ds
        logger.warning(f"tile {row},{col}: filling {layers} with nodata ({reason})")
        for name in layers:
            ds[name] = xr.full_like(template, np.nan)

    def _fetch(fn, *args, layers, coverage=None, **kwargs):
        nonlocal ds
        if coverage is not None and (tile_south > coverage[1] or tile_north < coverage[0]):
            _fill(layers, f"tile outside {fn.__name__} coverage {coverage}")
            return
        straddles = coverage is not None and (tile_south < coverage[0] or tile_north > coverage[1])
        t0 = time.time()
        try:
            ds = _with_transient_retries(fn, *args, log=log, **kwargs)
        except NoDataInBounds as e:  # the source's own definitive empty signal
            _fill(layers, f"{fn.__name__}: {e}")
        except Exception as e:
            if straddles:  # partially past the coverage edge: expected
                _fill(layers, f"{fn.__name__} at coverage edge: {e}")
            else:          # inside documented coverage: a REAL failure
                raise
        log(f"tile {row},{col}:   {', '.join(layers)}: {step_names.get(fn.__name__, fn.__name__)} "
            f"({time.time() - t0:.0f}s)")

    _fetch(add_dem, tile, ds, layers=('dem',))
    _fetch(add_chili, tile, ds, layers=('chili',), coverage=LAYER_LAT_COVERAGE['chili'])
    _fetch(add_snow_class, tile, ds, mask_nodata=mask_nodata, layers=('snow_classification',))
    _fetch(add_esa_worldcover, tile, ds, mask_nodata=mask_nodata, layers=('esa_worldcover',),
           coverage=LAYER_LAT_COVERAGE['esa_worldcover'])
    _fetch(add_forest_cover, tile, ds, mask_nodata=mask_nodata, layers=('forest_cover_fraction',),
           coverage=LAYER_LAT_COVERAGE['forest_cover_fraction'])
    # vector IDs: internal empty-handling only, exceptions must propagate
    t0 = time.time()
    ds = add_mountain_range_and_basin_and_continent(tile, ds, log=log)
    log(f"tile {row},{col}:   GMBA_V2_ID, PFAF_ID, continent: polygons rasterized (geocube 0.0003 deg) -> mode "
        f"({time.time() - t0:.0f}s)")
    ds = ds.drop_vars('_template')[list(ANCILLARY_LAYERS)]
    for name in ds.data_vars:   # keep only serializable attrs (easysnowdata >= 0.0.25 attaches dicts)
        ds[name].attrs = {k: v for k, v in ds[name].attrs.items()
                          if isinstance(v, (str, bytes, int, float, list, tuple, np.ndarray, np.number))}
    ds.attrs.update({'tile_row': row, 'tile_col': col, 'grid': f"{config.version}_grid",
                     'basin_atlas_layer': settings.BASIN_ATLAS_LAYER})
    return ds


def _encode_layer(name, values):
    """Decoded layer values -> the stored integer array (CF inverse: (v - offset) / scale,
    nodata -> the encoded fill)."""
    e = ANCILLARY_ENCODING[name]
    v = np.asarray(values, dtype='float64')
    v = np.where(v == -9999, np.nan, v)                  # the builders' integer nodata sentinel
    if 'scale_factor' in e:
        v = (v - e.get('add_offset', 0.0)) / e['scale_factor']
    info = np.iinfo(e['dtype'])
    v = np.where((v < info.min) | (v > info.max), np.nan, v)   # never let a value wrap (int8 -9999 -> -15, 2026-09-04)
    out = np.where(np.isfinite(v), np.round(v), e['_FillValue'])
    return out.astype(e['dtype'])


def write_ancillary_tile(config, repo, row, col, ds, layers=ANCILLARY_LAYERS, duration_s=None, log=None):
    """Write a tile's layers into their chunks of the ancillary store and commit them as ONE
    ledger entry (kind ``ancillary_tile``; newest wins, so a rewrite or a unit-layer refresh
    supersedes the previous commit). Returns the snapshot id."""
    import zarr
    from global_snowmelt_runoff_onset import store as gs_store
    log = log or logger.info
    region = gs_store.tile_region_slices(config, row, col)
    ys, xs = region['latitude'], region['longitude']
    shape = (ys.stop - ys.start, xs.stop - xs.start)
    encoded, stats = {}, {}
    for name in layers:
        values = ds[name].values
        if values.shape != shape:
            raise ValueError(f"tile {row},{col}: {name} has shape {values.shape}, the tile window is {shape}")
        encoded[name] = _encode_layer(name, values)
        stats[name] = round(float(np.mean(encoded[name] != ANCILLARY_ENCODING[name]['_FillValue'])), 4)

    def write_fn(session):
        g = zarr.open_group(session.store, mode='r+')
        for name, arr in encoded.items():
            g[name][ys, xs] = arr

    metadata = {'schema': ledger.SCHEMA, 'kind': KIND_ANCILLARY_TILE, 'tile': [int(row), int(col)],
                'grid': f"{config.version}_grid", 'layers': list(layers), 'valid_fraction': stats,
                'basin_atlas_layer': settings.BASIN_ATLAS_LAYER,
                'duration_s': round(float(duration_s), 1) if duration_s is not None else None,
                'provenance': ledger.provenance()}
    snap = ledger.commit_with_retry(repo, write_fn, f"tile {row:03d}_{col:03d}: ancillary ({len(layers)} layers)",
                                    metadata, log=log)
    log(f"tile {row},{col}:   {len(layers)} layers written to the ancillary store and committed -> {snap}")
    return snap


def open_ancillary_window(config, row, col, repo=None, local_store=None):
    """The tile's window of the ancillary store, decoded (nodata -> NaN, CHILI in 0-1) and
    georeferenced (EPSG:4326). An uncommitted tile reads back as all-NaN, so callers check
    the ledger first."""
    from global_snowmelt_runoff_onset import store as gs_store
    repo = repo or open_ancillary_repo(config, local_store)
    ds = xr.open_zarr(repo.readonly_session(ledger.BRANCH).store, consolidated=False, zarr_format=3,
                      decode_coords='all', chunks=None)
    window = ds.isel(gs_store.tile_region_slices(config, row, col)).load()
    # a continent code below 0 is nodata: tiles built before 2026-09-04 stored the -9999 "no polygon"
    # fill wrapped to -15 in int8 (see _encode_layer)
    window['continent'] = window['continent'].where(window['continent'] >= 0)
    return window.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude').rio.write_crs('EPSG:4326')


def refresh_unit_layers(config, row, col, layers=('GMBA_V2_ID', 'PFAF_ID', 'continent'), repo=None, log=None):
    """Stage 0b: re-rasterize the vector unit-id layers of a STORED tile into the ancillary
    store (one commit; the raster layers and their Earth Engine work are untouched). How a
    unit definition changes after the ancillary exists (e.g. another HydroBASINS level);
    the tile's partials must be re-mapped afterwards (:func:`process_tile`)."""
    log = log or logger.info
    repo = repo or open_ancillary_repo(config)
    tile = config.get_tile(row, col)
    window = open_ancillary_window(config, row, col, repo)
    ds = window[['dem']].copy()
    ds = add_mountain_range_and_basin_and_continent(tile, ds, log=log)
    snap = write_ancillary_tile(config, repo, row, col, ds, layers=layers, log=log)
    check = open_ancillary_window(config, row, col, repo)
    for name in layers:
        if not np.array_equal(_encode_layer(name, check[name].values), _encode_layer(name, ds[name].values)):
            raise RuntimeError(f"tile {row},{col}: {name} did not read back as written")
    log(f"tile {row},{col}: refreshed {layers} ({snap})")
    return check[list(layers)]


def fresh_blob_fs(config):
    """See :func:`gsro_analysis.settings.fresh_blob_fs`."""
    return settings.fresh_blob_fs(config)


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


# ---------------------------------------------------------------------------
# Process tiles to parquets: both stores' windows -> the 80 m UTM grid -> pixel table -> partial sums

def terrain_derivatives(dem_utm):
    """Aspect and slope (degrees) of a DEM on a metric grid (xdem, resolution from the grid)."""
    attributes = xdem.terrain.get_terrain_attribute(dem_utm, resolution=abs(float(dem_utm.rio.resolution()[0])),
                                                    attribute=["aspect", "slope"])
    aspect = xr.DataArray(attributes[0], dims=dem_utm.dims, coords=dem_utm.coords)
    slope = xr.DataArray(attributes[1], dims=dem_utm.dims, coords=dem_utm.coords)
    return aspect, slope


def tabulate_tile(config, row, col, ancillary_window, global_ds=None):
    """Join one tile's runoff-onset values (any dataset version) and its ancillary window on
    the tile's 80 m UTM grid and tabulate.

    Both windows come off the dataset grid (the same pixel indices); each is reprojected onto
    :func:`tile_utm_template` (onset, DEM, CHILI, forest cover bilinear; the categorical and
    id layers nearest, :data:`UTM_RESAMPLING`), then slope and aspect are derived from the
    UTM DEM. Rows are UTM pixels (~6,400 m2), so counts are area. The store is read with
    mask_and_scale=True and chunks=None (per-tile reads are small; CF decoding rescales
    runoff_onset_mad / temporal_resolution_median and turns -9999 into NaN before the
    bilinear reprojection). An xr.merge(join='exact') tripwire guards the alignment.
    """
    from global_snowmelt_runoff_onset import store as gs_store

    if global_ds is None:
        global_ds = config.open_runoff_onset_dataset(chunks=None, mask_and_scale=True)
    template = tile_utm_template(config, row, col)
    region = gs_store.tile_region_slices(config, row, col)
    tile_ds = global_ds.isel(region).drop_vars('temporal_resolution')
    tile_ds = tile_ds.astype(np.float32).compute()

    onset_utm = tile_ds.rio.reproject_match(template, resampling=rasterio.enums.Resampling.bilinear)
    onset_utm = convert_water_year_dim_to_var(onset_utm)

    anc = xr.Dataset()
    for name in ANCILLARY_LAYERS:
        da = ancillary_window[name].astype(np.float32)
        if name in ('GMBA_V2_ID', 'PFAF_ID'):
            da = da.fillna(-9999.0)          # nodata as the integer sentinel tabulate expects
        anc[name] = da.rio.write_nodata(np.nan).rio.reproject_match(template, resampling=UTM_RESAMPLING[name])
    anc['aspect'], anc['slope'] = terrain_derivatives(anc['dem'])
    # exact WGS84 pixel-center coordinates of the UTM grid
    anc['original_lat'], anc['original_lon'] = _utm_latlon_arrays(anc['dem'])
    joined = xr.merge([onset_utm, anc], join='exact', compat='no_conflicts',
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
    df.attrs['utm_crs'] = str(template.rio.crs)
    df.attrs['utm_shape'] = (int(template.sizes['y']), int(template.sizes['x']))
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


def tabulate_and_save_tile(config, row, col, global_ds=None, repo=None):
    """Tabulate one tile and write its pixel table (opt-in product)."""
    window = open_ancillary_window(config, row, col, repo)
    df = tabulate_tile(config, row, col, window, global_ds=global_ds)
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


def write_partials(df, config, row, col, filter_tags=None, log=None):
    log = log or logger.info
    """Partial sums of a tile's pixel table -> partials_tile_path() (single
    blob). Returns the number of partial rows."""
    from gsro_analysis import aggregate

    partials = aggregate.tile_partials(df, config.water_years, filter_tags=filter_tags)
    partials.insert(0, 'tile_col', np.int16(col))
    partials.insert(0, 'tile_row', np.int16(row))
    path = partials_tile_path(config, row, col)
    partials.to_parquet(path, filesystem=config.azure_blob_fs, index=False,
                        compression='zstd')
    if len(partials):
        log(f"tile {row},{col}: map — {len(partials):,} partial rows "
            f"({partials['filter_tag'].nunique()} filters x {partials['unit_type'].nunique()} unit types) -> {path}")
    else:   # a verified-empty tile (no pixel with a valid median passed the filters): the empty blob IS the result
        log(f"tile {row},{col}: map — 0 partial rows (no pixel with a valid median) -> empty {path}")
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
    'ancillary_grid': lambda c: f"{settings.ANCILLARY_PREFIX}/{c.version}_grid",    # the ancillary icechunk repo (per GRID; rebuilding = Earth Engine again)
    'aggregated_mirror': lambda c: f"{settings.AGGREGATED_PREFIX}/{c.version}",    # reduce --mirror copies of the cubes
    'era5_land':     lambda c: f"{settings.ERA5_LAND_PREFIX}/{c.version}",          # the ERA5-Land icechunk repo (acquisition + anomaly group;
                                                                                    # Get ERA5-Land data rebuilds it, ~1 h) — NOT in the default reset
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
    'Get ERA5-Land data' workflow re-acquires every water year; its start_fresh box does the same). Never deletes the icechunk dataset,
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
    print("done; the next 'Get ancillary data' / 'Process tiles to parquets' dispatch rebuilds every tile")
    return plan


def process_tile(config, row, col, global_ds=None, keep_pixels=False,
                 filter_tags=None, log=None, repo=None):
    """The whole per-tile job of 'Process tiles to parquets': the ancillary window (must have
    a commit) and the store window -> the 80 m UTM grid -> pixel table -> partial sums
    written to partials_tile_path(); the pixel table itself only with ``keep_pixels``.
    Failure = exception = no partials blob (the fleet rule). ``log`` (default: the module
    logger) receives one line per step. Returns (n_pixels, n_partial_rows)."""
    log = log or logger.info
    repo = repo or open_ancillary_repo(config)
    if not ancillary_tile_complete(config, row, col, repo):
        raise RuntimeError(f"tile {row},{col}: no ancillary commit in {ancillary_repo_prefix(config)}: "
                           "run the 'Get ancillary data' workflow first")
    t0 = time.time()
    window = open_ancillary_window(config, row, col, repo)
    log(f"tile {row},{col}: ancillary window read from the dataset-grid store "
        f"({window.sizes['latitude']} x {window.sizes['longitude']} px, {len(ANCILLARY_LAYERS)} layers) "
        f"({time.time() - t0:.0f}s)")
    t1 = time.time()
    df = tabulate_tile(config, row, col, window, global_ds=global_ds)
    shape = df.attrs.get('utm_shape', ('?', '?'))
    log(f"tile {row},{col}: tabulate: store window ({len(config.water_years)} water years + medians) and "
        f"ancillary window reprojected onto the {df.attrs.get('utm_crs')} 80 m grid "
        f"({shape[0]} x {shape[1]} px), slope + aspect derived, {len(df):,} pixels with a valid median "
        f"({time.time() - t1:.0f}s)")
    if keep_pixels:
        save_pixel_table(df, config, row, col)
        log(f"tile {row},{col}:   pixel table written -> {parquet_tile_path(config, row, col)}")
    n_rows = write_partials(df, config, row, col, filter_tags=filter_tags, log=log)
    return len(df), n_rows
