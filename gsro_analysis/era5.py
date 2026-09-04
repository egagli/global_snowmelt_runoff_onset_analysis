"""ERA5-Land for the analyses: ONE icechunk repository per dataset version.

Layout on Azure — ``settings.ERA5_LAND_PREFIX/<version>/era5_land`` is an icechunk
repository (Zarr v3 inside; open it with :func:`open_era5_land`, never as a plain
zarr path):

    /            the acquisition: the 8 ERA5-Land monthly variables (:data:`VARIABLES`) on
                 ``(water_year, month, latitude, longitude)`` — every water year of
                 ``config.water_years``, 9 hemisphere-aware months (:data:`MONTH_LABELS`),
                 the NATIVE 0.1° grid: 1800 x 3600 cells centred at multiples of 0.1°,
                 north-down (:data:`LATITUDE` 89.9..-90.0, :data:`LONGITUDE` -180.0..179.9;
                 the +90° row is Arctic Ocean, which ERA5-Land does not cover, and the
                 +180° column duplicates -180°)
    /anomaly     the same layout minus the per-pixel median over ALL water years (the
                 base period is in the group attrs and in the commit metadata)

The repository IS the ledger (the production repo's icechunk fleet pattern,
``docs/icechunk-github-actions-pattern.md`` there): the store is initialized
once as an empty template (:func:`initialize`), every water year is ONE commit
carrying machine-readable metadata (:func:`write_water_year`; one chunk = one
variable x water year x month, so concurrent year jobs never touch the same
chunk and their commits rebase automatically), the anomaly is one commit
(:func:`build_anomaly`), and "what remains" is a fold over the commit history
(:func:`status`). A failed job commits nothing and its year is simply re-listed
by the next dispatch. Acquisitions are never copied between dataset versions.

Acquisition goes through xee >= 0.1 with the collection's native grid given
explicitly (:data:`NATIVE_GRID`), so every value is the asset's own, never a
resample: the pre-2026-09-03 per-year stores were a nearest-neighbour resample
half a cell off this grid. One variable x hemisphere window is in memory at a
time (~230 MB); a water year takes about a minute.

Needs Earth Engine (``settings.initialize_earthengine()``, service-account key,
works headless) and the Azure SAS token (through the production ``Config``).
"""

import time
import warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr

from gsro_analysis import ledger, settings

EE_COLLECTION = "ECMWF/ERA5_LAND/MONTHLY_AGGR"

# 9 hemisphere-aware months per water year: NH Dec(wy-1)-Aug(wy), SH Jun(wy)-Feb(wy+1),
# both labeled winter -> spring -> summer; cells with latitude >= 0 carry the NH window.
MONTH_LABELS = ['winter_month_1', 'winter_month_2', 'winter_month_3',
                'spring_month_1', 'spring_month_2', 'spring_month_3',
                'summer_month_1', 'summer_month_2', 'summer_month_3']

VARIABLES = ["temperature_2m",
             "dewpoint_temperature_2m",
             "snow_depth_water_equivalent",
             "snowfall_sum",
             "surface_latent_heat_flux_sum",
             "surface_sensible_heat_flux_sum",
             "surface_solar_radiation_downwards_sum",
             "surface_thermal_radiation_downwards_sum"]

# units per the ERA5-Land monthly aggregates catalog entry (the bands carry none)
VARIABLE_ATTRS = {
    "temperature_2m": {"units": "K", "long_name": "2 m air temperature, monthly mean"},
    "dewpoint_temperature_2m": {"units": "K", "long_name": "2 m dewpoint temperature, monthly mean"},
    "snow_depth_water_equivalent": {"units": "m of water equivalent",
                                    "long_name": "snow depth water equivalent, monthly mean"},
    "snowfall_sum": {"units": "m of water equivalent", "long_name": "snowfall, monthly sum of daily totals"},
    "surface_latent_heat_flux_sum": {"units": "J m-2", "long_name": "surface latent heat flux, monthly sum"},
    "surface_sensible_heat_flux_sum": {"units": "J m-2", "long_name": "surface sensible heat flux, monthly sum"},
    "surface_solar_radiation_downwards_sum": {"units": "J m-2",
                                              "long_name": "surface solar radiation downwards, monthly sum"},
    "surface_thermal_radiation_downwards_sum": {"units": "J m-2",
                                                "long_name": "surface thermal radiation downwards, monthly sum"},
}
DATA_CITATION = ("Muñoz Sabater, J. (2019): ERA5-Land monthly averaged data from 1950 to present. "
                 "Copernicus Climate Change Service (C3S) Climate Data Store (CDS), doi:10.24381/cds.68d2bb30; "
                 "accessed through Google Earth Engine " + EE_COLLECTION)

# the native grid (Earth Engine's crs_transform for the collection is
# [0.1, 0, -180.05, 0, -0.1, 90.05] over 3601 x 1801 pixels; see the module docstring)
RES = 0.1
N_LAT, N_LON = 1800, 3600
LATITUDE = np.round(89.9 - RES * np.arange(N_LAT), 1)       # north-down, cell centres
LONGITUDE = np.round(-180.0 + RES * np.arange(N_LON), 1)
CRS_TRANSFORM = [RES, 0, -180.05, 0, -RES, 89.95]           # upper-left corner of the stored grid
DIMS = ("water_year", "month", "latitude", "longitude")
CHUNKS = (1, 1, N_LAT, N_LON)                               # one chunk = one variable x year x month (~26 MB raw)

ANOMALY_GROUP = "anomaly"
BRANCH = "main"
SCHEMA = 1
KIND_INIT, KIND_WATER_YEAR, KIND_ANOMALY = "init", "water_year", "anomaly"
COMMIT_MAX_TRIES = 6


# the collection's native grid as xee >= 0.1 wants it (== xee.helpers.extract_grid_params(collection)):
# 3601 x 1801 pixels, row 0 = latitude 90.0 (dropped: Arctic Ocean), column 3600 = +180.0 (dropped: duplicate)
NATIVE_GRID = {'crs': 'EPSG:4326', 'crs_transform': (RES, 0, -180.05, 0, -RES, 90.05), 'shape_2d': (N_LON + 1, N_LAT + 1)}


# ---------------------------------------------------------------------------
# the repository

def repo_prefix(config):
    """Container-qualified Azure prefix of the version's icechunk repository."""
    return f"{settings.ERA5_LAND_PREFIX}/{config.version}/era5_land"


def repo_storage(config):
    return ledger.azure_storage(config, repo_prefix(config))


def repo_config():
    return ledger.repo_config()


def _storage(config, local_store=None):
    return ledger.storage(config, repo_prefix(config), local_store)


def repo_exists(config, local_store=None):
    return ledger.repo_exists(config, repo_prefix(config), local_store)


def open_repo(config, local_store=None):
    """Open the version's repository (``local_store``: a local icechunk repo path, for tests);
    a clear FileNotFoundError when it has not been created yet."""
    return ledger.open_repo(config, repo_prefix(config), local_store, what='ERA5-Land repository')


def _provenance():
    return ledger.provenance()


def build_template(config):
    """The lazy, all-NaN template of the acquisition (and of the anomaly group) plus its
    Zarr v3 encoding; written metadata-only by :func:`initialize`."""
    import dask.array as da
    import rioxarray  # noqa: F401
    import zarr

    years = [int(y) for y in config.water_years]
    shape = (len(years), len(MONTH_LABELS), N_LAT, N_LON)
    ds = xr.Dataset(
        {v: (DIMS, da.full(shape, np.nan, dtype='float32', chunks=CHUNKS), dict(VARIABLE_ATTRS[v]))
         for v in VARIABLES},
        coords={'water_year': years, 'month': MONTH_LABELS, 'latitude': LATITUDE, 'longitude': LONGITUDE})
    ds['latitude'].attrs = {'units': 'degrees_north', 'long_name': 'latitude of the cell centre (native ERA5-Land 0.1 degree grid, north-down)'}
    ds['longitude'].attrs = {'units': 'degrees_east', 'long_name': 'longitude of the cell centre (native ERA5-Land 0.1 degree grid)'}
    ds['month'].attrs = {'description': ('hemisphere-aware month of the water year: cells with latitude >= 0 carry '
                                         'Dec(wy-1)..Aug(wy), cells south of the equator Jun(wy)..Feb(wy+1)')}
    ds['water_year'].attrs = {'description': 'water year (NH Oct(wy-1)-Sep(wy); SH Apr(wy)-Mar(wy+1))'}
    ds = ds.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude').rio.write_crs('EPSG:4326')
    for v in VARIABLES:                     # an explicit attr: rioxarray keeps grid_mapping in encoding, which the
        ds[v].attrs['grid_mapping'] = 'spatial_ref'   # zarr write drops, and decode_coords='all' needs it
    ds.attrs.update({
        'title': 'ERA5-Land monthly aggregates for the global snowmelt runoff onset analyses',
        'source': f'Google Earth Engine {EE_COLLECTION} through xee on the native grid (no resampling)',
        'crs': 'EPSG:4326', 'crs_transform': CRS_TRANSFORM, 'cadence': 'MONTHLY',
        'data_citation': DATA_CITATION, 'dataset_version': config.version,
    })
    compressor = zarr.codecs.BloscCodec(cname='zstd', clevel=5)
    encoding = {v: {'chunks': CHUNKS, 'compressors': [compressor], 'dtype': 'float32',
                    '_FillValue': np.nan, 'fill_value': np.nan} for v in VARIABLES}
    return ds, encoding


def initialize(config, start_fresh=False, local_store=None, log=print):
    """Return the version's repository, creating it with the empty template if it does not
    exist. ``start_fresh=True`` DELETES an existing repository first (the 'Get ERA5-Land data'
    workflow's off-by-default box; the only deletion in this module)."""
    if repo_exists(config, local_store):
        if not start_fresh:
            log(f"repository exists: {local_store or repo_prefix(config)}")
            return open_repo(config, local_store)
        ledger.delete_repo(config, repo_prefix(config), local_store, log=log)
    years = [int(y) for y in config.water_years]
    repo = ledger.create_repo(config, repo_prefix(config), local_store)
    template, encoding = build_template(config)
    session = repo.writable_session(BRANCH)
    template.to_zarr(session.store, mode='w', zarr_format=3, compute=False, write_empty_chunks=False,
                     consolidated=False, encoding=encoding)
    session.commit(f"initialize empty ERA5-Land store, WY{years[0]}-{years[-1]}",
                   metadata={'schema': SCHEMA, 'kind': KIND_INIT, 'config_version': config.version,
                             'water_years': years, 'variables': list(VARIABLES),
                             'grid': {'n_lat': N_LAT, 'n_lon': N_LON, 'crs_transform': CRS_TRANSFORM},
                             'provenance': _provenance()})
    log(f"initialized {local_store or repo_prefix(config)}: empty template, WY{years[0]}-{years[-1]}")
    return repo


# ---------------------------------------------------------------------------
# the ledger: a fold over the commit history

def commit_records(repo, branch=BRANCH):
    return ledger.commit_records(repo, branch)


def status(config, repo=None, local_store=None):
    """What the repository holds for ``config.water_years``: ``years`` {wy: newest commit
    record}, ``remaining`` (missing water years, in order), ``complete``, ``anomaly`` (its
    commit record or None), ``anomaly_stale`` (a year committed after it, a year missing,
    or a different base period) and the branch tip ``snapshot_id``."""
    repo = repo or open_repo(config, local_store)
    expected = [int(y) for y in config.water_years]
    years, anomaly = {}, None
    for r in commit_records(repo):                       # newest first: first seen wins
        if r['kind'] == KIND_WATER_YEAR:
            years.setdefault(int(r['water_year']), r)
        elif r['kind'] == KIND_ANOMALY and anomaly is None:
            anomaly = r
    remaining = [y for y in expected if y not in years]
    newest_year = min((r['ancestry_index'] for y, r in years.items() if y in expected), default=None)
    stale = anomaly is not None and (
        bool(remaining) or (newest_year is not None and newest_year < anomaly['ancestry_index'])
        or [int(y) for y in anomaly.get('base_water_years', [])] != expected)
    return {'years': years, 'remaining': remaining, 'complete': not remaining,
            'anomaly': anomaly, 'anomaly_stale': stale, 'snapshot_id': repo.lookup_branch(BRANCH)}


def commit_with_retry(repo, write_fn, message, metadata, branch=BRANCH, allow_empty=False,
                      max_tries=COMMIT_MAX_TRIES, log=print):
    return ledger.commit_with_retry(repo, write_fn, message, metadata, branch=branch, allow_empty=allow_empty,
                                    max_tries=max_tries, log=log)


# ---------------------------------------------------------------------------
# acquisition (Earth Engine, native grid)

def month_windows(water_year):
    """[(label, (NH year, month), (SH year, month))] for the 9 months of a water year."""
    nh = [(water_year - 1, 12)] + [(water_year, m) for m in range(1, 9)]
    sh = [(water_year, m) for m in range(6, 13)] + [(water_year + 1, 1), (water_year + 1, 2)]
    return list(zip(MONTH_LABELS, nh, sh))


def _open_window(variables, start, end):
    """The collection's monthly images in [start, end) as a lazy xee Dataset (time, y, x)
    on the native grid, time ascending. Earth Engine must be initialized."""
    import ee
    import xee  # noqa: F401  -- registers engine='ee'
    ic = ee.ImageCollection(EE_COLLECTION).filterDate(start, end).select(list(variables))
    return xr.open_dataset(ic, engine='ee', **NATIVE_GRID).sortby('time')


def fetch_water_year(water_year, variables=VARIABLES, log=print):
    """{variable: float32 (9, 1800, 3600)} for one water year from Earth Engine through xee
    on the native grid: the NH window (Dec(wy-1)..Aug(wy)) fills the rows with latitude
    >= 0, the SH window (Jun(wy)..Feb(wy+1)) the rows south of the equator; NaN where
    ERA5-Land has no data. Raises if a window does not hold its 9 monthly images yet
    (failure = nothing written). One variable x hemisphere (~230 MB) in memory at a time."""
    wy = int(water_year)
    windows = {'NH': (f"{wy - 1}-12-01", f"{wy}-09-01"), 'SH': (f"{wy}-06-01", f"{wy + 1}-03-01")}
    expected_months = {'NH': [12, 1, 2, 3, 4, 5, 6, 7, 8], 'SH': [6, 7, 8, 9, 10, 11, 12, 1, 2]}
    half = N_LAT // 2
    out = {v: np.full((len(MONTH_LABELS), N_LAT, N_LON), np.nan, dtype='float32') for v in variables}
    t0 = time.time()
    for hemi, (start, end) in windows.items():
        ds = _open_window(variables, start, end)
        months = [int(m) for m in ds['time'].dt.month.values]
        if months != expected_months[hemi]:
            raise RuntimeError(f"{EE_COLLECTION}: the {hemi} window [{start}, {end}) of WY{wy} holds months "
                               f"{months}, expected {expected_months[hemi]} (not all months published yet?)")
        rows = slice(1, 1 + half) if hemi == 'NH' else slice(1 + half, 1 + N_LAT)   # skip the +90 row
        for v in variables:
            arr = ds[v].isel(y=rows, x=slice(0, N_LON)).values                    # (9, 900, 3600), float64
            out[v][:, rows.start - 1:rows.stop - 1] = np.where(np.isfinite(arr), arr, np.nan).astype('float32')
            del arr
        log(f"  WY{wy} {hemi} window: {len(variables)} variables ({time.time() - t0:.0f}s)")
    return out


def _qa_stats(arrays):
    """Small QA numbers for the commit metadata: finite fraction and land mean of the
    first spring month of every variable."""
    mi = MONTH_LABELS.index('spring_month_1')
    stats = {}
    for v, arr in arrays.items():
        slab = arr[mi]
        finite = np.isfinite(slab)
        stats[v] = {'finite_fraction': round(float(finite.mean()), 4),
                    'mean': round(float(slab[finite].mean()), 4) if finite.any() else None}
    return stats


def write_water_year(config, repo, water_year, arrays, duration_s=None, log=print):
    """Write one water year's arrays into their slots and commit it as ONE ledger entry
    (kind ``water_year``). Refuses arrays with no land data. Superseding an existing
    year's commit is allowed (newest wins) but logged loudly."""
    import zarr
    water_year = int(water_year)
    years = [int(y) for y in config.water_years]
    j = years.index(water_year)
    stats = _qa_stats(arrays)
    frac = stats['temperature_2m']['finite_fraction'] if 'temperature_2m' in stats else max(s['finite_fraction'] for s in stats.values())
    if frac < 0.1:
        raise RuntimeError(f"WY{water_year}: only {frac:.3f} of the spring slab is finite — refusing to write")
    if water_year in status(config, repo)['years']:
        log(f"WARNING: WY{water_year} already has a commit; this write supersedes it")

    def write_fn(session):
        g = zarr.open_group(session.store, mode='r+')
        assert int(g['water_year'][j]) == water_year, "water_year coordinate does not match the store"
        for v, arr in arrays.items():
            g[v][j] = arr

    metadata = {'schema': SCHEMA, 'kind': KIND_WATER_YEAR, 'water_year': water_year,
                'config_version': config.version, 'variables': sorted(arrays), 'stats': stats,
                'duration_s': round(float(duration_s), 1) if duration_s is not None else None,
                'provenance': _provenance()}
    snap = commit_with_retry(repo, write_fn, f"WY{water_year}: ERA5-Land monthly, {len(arrays)} variables",
                             metadata, log=log)
    log(f"committed WY{water_year} -> {snap}")
    return snap


def acquire_water_year(config, water_year, repo=None, local_store=None, log=print):
    """Fetch + write + commit one water year (the fleet job's whole work)."""
    repo = repo or open_repo(config, local_store)
    t0 = time.time()
    arrays = fetch_water_year(water_year, log=log)
    return write_water_year(config, repo, water_year, arrays, duration_s=time.time() - t0, log=log)


# ---------------------------------------------------------------------------
# opening + the anomaly group

def open_era5_land(config, repo=None, local_store=None, chunks='auto', group=None):
    """The acquisition (or, with ``group='anomaly'``, the anomaly) as a lazy Dataset on
    ``(water_year, month, latitude, longitude)``, north-down native grid, CF coords."""
    repo = repo or open_repo(config, local_store)
    ds = xr.open_zarr(repo.readonly_session(BRANCH).store, group=group, zarr_format=3,
                      consolidated=False, decode_coords='all', chunks=chunks)
    if 'spatial_ref' in ds.data_vars:      # stores written before grid_mapping was an explicit attr
        ds = ds.set_coords('spatial_ref')
    return ds


def build_anomaly(config, repo=None, local_store=None, force=False, log=print):
    """(Re)build the ``anomaly`` group — every variable minus its per-pixel median over
    all water years — once every water year is committed; skipped while years are
    missing, and when an up-to-date anomaly commit exists (unless ``force``). One
    (variable, month) slab at a time (~1 GB peak), one commit (kind ``anomaly``)."""
    import zarr
    repo = repo or open_repo(config, local_store)
    st = status(config, repo)
    if st['remaining']:
        log(f"anomaly not built: water years still missing {st['remaining']}")
        return None
    if st['anomaly'] is not None and not st['anomaly_stale'] and not force:
        log(f"anomaly up to date (commit {st['anomaly']['snapshot_id']})")
        return st['anomaly']['snapshot_id']
    years = [int(y) for y in config.water_years]
    src = open_era5_land(config, repo, chunks=None)
    template, encoding = build_template(config)
    template.attrs.update({'title': 'ERA5-Land monthly anomalies vs the per-pixel median over the base water years',
                           'anomaly_base_water_years': years})
    for v in VARIABLES:
        template[v].attrs['long_name'] += ', anomaly vs the median over all water years'
    t0 = time.time()

    def write_fn(session):
        template.to_zarr(session.store, group=ANOMALY_GROUP, mode='w', zarr_format=3, compute=False,
                         write_empty_chunks=False, consolidated=False, encoding=encoding)
        g = zarr.open_group(session.store, mode='r+')[ANOMALY_GROUP]
        for v in VARIABLES:
            for mi in range(len(MONTH_LABELS)):
                slab = src[v].isel(month=mi).values                    # (n_years, 1800, 3600)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)      # all-NaN ocean columns
                    median = np.nanmedian(slab, axis=0)
                g[v][:, mi] = (slab - median).astype('float32')
                del slab
            log(f"  anomaly {v} done ({time.time() - t0:.0f}s)")

    metadata = {'schema': SCHEMA, 'kind': KIND_ANOMALY, 'config_version': config.version,
                'base_water_years': years, 'source_snapshot_id': st['snapshot_id'],
                'variables': list(VARIABLES), 'provenance': _provenance()}
    snap = commit_with_retry(repo, write_fn, f"anomaly vs the median over WY{years[0]}-{years[-1]}",
                             metadata, max_tries=3, log=log)
    log(f"committed anomaly -> {snap} ({time.time() - t0:.0f}s)")
    return snap


def open_anomaly(config, repo=None, local_store=None, chunks='auto'):
    """The monthly anomaly group feeding the per-unit temperature sensitivity; refuses a
    repository without an anomaly commit and warns when it is stale."""
    repo = repo or open_repo(config, local_store)
    st = status(config, repo)
    if st['anomaly'] is None:
        raise FileNotFoundError(
            f"no anomaly commit in {repo_prefix(config)}: run the 'Get ERA5-Land data' workflow (it builds the "
            "anomaly once every water year is committed) or era5.build_anomaly(config)")
    if st['anomaly_stale']:
        warnings.warn("the ERA5-Land anomaly group is stale (a water year was committed after it, or the "
                      "base period differs from config.water_years): rebuild it", stacklevel=2)
    return open_era5_land(config, repo, chunks=chunks, group=ANOMALY_GROUP)


def correlations_store_path(config):
    """LOCAL output of the exploratory pixel-wise correlation notebook
    (analyses/climate/pixelwise_climate_correlations.ipynb): a zarr under the
    gitignored scratch/ folder, keyed by dataset version. Its own product,
    read by nothing else, so it never goes to Azure (2026-09-03)."""
    from gsro_analysis import paths
    paths.SCRATCH.mkdir(parents=True, exist_ok=True)
    return paths.SCRATCH / f"{config.version}_runoff_onset_and_era5_eqearth_anomaly_correlations.zarr"


# ---------------------------------------------------------------------------
# zonal anomalies per analysis unit (replaces the per-range rio.clip join
# that used to live in stage 2 and needed a 16 GB machine)

def _era5_grid(anomaly_ds):
    """(lat ascending?, lat, lon, cell size) of the native 0.1 deg store."""
    lat = anomaly_ds["latitude"].values
    lon = anomaly_ds["longitude"].values
    return lat, lon, float(abs(lat[1] - lat[0]))


def unit_coverage_matrix(units_gdf, id_col, anomaly_ds, subcells=5,
                         geometry_fixes=None, skip_names=()):
    """Sparse (n_units x n_cells) matrix of the FRACTION of each ERA5 cell
    covered by each unit polygon, from a ``subcells`` x ``subcells``
    sub-sampling of every cell (5 -> 0.02 deg, ~2 km). Cells are flattened
    in the store's (latitude, longitude) order. Returns (matrix, unit_ids).
    ``geometry_fixes``: {name: bbox} pre-clips (aggregate.RANGE_GEOMETRY_FIXES)."""
    import rasterio.features
    import rasterio.transform
    import scipy.sparse as sp

    lat, lon, res = _era5_grid(anomaly_ds)
    nlat, nlon = len(lat), len(lon)
    lat_asc = lat[1] > lat[0]
    lat_top = (lat.max() + res / 2)
    lon_left = (lon.min() - res / 2)
    sub = res / subcells
    rows, cols, vals, ids = [], [], [], []
    name_col = "MapName" if "MapName" in units_gdf else None
    for k, (_, unit) in enumerate(units_gdf.iterrows()):
        name = unit[name_col] if name_col else None
        if name in skip_names:
            continue
        geom = unit.geometry
        if geometry_fixes and name in geometry_fixes:
            from shapely.geometry import box
            geom = geom.intersection(box(*geometry_fixes[name]))
        if geom.is_empty:
            continue
        minx, miny, maxx, maxy = geom.bounds
        # window of whole ERA5 cells (north-up rows) containing the bounds
        c0 = int(np.floor((minx - lon_left) / res)); c1 = int(np.ceil((maxx - lon_left) / res))
        r0 = int(np.floor((lat_top - maxy) / res)); r1 = int(np.ceil((lat_top - miny) / res))
        c0, c1 = max(c0, 0), min(c1, nlon); r0, r1 = max(r0, 0), min(r1, nlat)
        if c1 <= c0 or r1 <= r0:
            continue
        transform = rasterio.transform.from_origin(lon_left + c0 * res, lat_top - r0 * res, sub, sub)
        shape = ((r1 - r0) * subcells, (c1 - c0) * subcells)
        mask = rasterio.features.rasterize([(geom, 1)], out_shape=shape, transform=transform,
                                           fill=0, dtype="uint8", all_touched=False)
        frac = mask.reshape(r1 - r0, subcells, c1 - c0, subcells).sum(axis=(1, 3)) / subcells ** 2
        rr, cc = np.nonzero(frac)
        if rr.size == 0:
            continue
        row_north_up = rr + r0
        lat_idx = (nlat - 1 - row_north_up) if lat_asc else row_north_up
        rows.append(np.full(rr.size, len(ids))); cols.append(lat_idx * nlon + cc + c0); vals.append(frac[rr, cc])
        ids.append(unit[id_col])
    W = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(ids), nlat * nlon))
    return W, np.array(ids)


def pyramid_masks(config, anomaly_ds, level=7):
    """Seasonal-snow fraction (lat, lon) and per-water-year onset validity
    (water_year, lat, lon) on the ERA5 grid, from the public pyramid level
    ``level`` (~10 km): ``seasonal_snow_pct`` is the Sturm & Liston
    non-ephemeral area fraction, validity = a runoff_onset value exists in
    the cell. These are the same two masks the old per-range clip applied
    (snow class != 4, onset notnull), now as ERA5-cell weights."""
    import rioxarray  # noqa: F401
    from global_snowmelt_runoff_onset.pyramid import open_pyramid_level

    pyr = open_pyramid_level(config, level)
    template = xr.DataArray(np.zeros((anomaly_ds.sizes["latitude"], anomaly_ds.sizes["longitude"]), dtype="float32"),
                            dims=("latitude", "longitude"),
                            coords={"latitude": anomaly_ds["latitude"], "longitude": anomaly_ds["longitude"]})
    template = template.rio.write_crs("EPSG:4326").rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
    import rasterio.enums as re_
    snow = pyr["seasonal_snow_pct"].rio.write_crs("EPSG:4326").astype("float32").fillna(0.0)
    snow_frac = (snow.rio.reproject_match(template, resampling=re_.Resampling.average) / 100.0).clip(0, 1)
    valid = []
    for wy in pyr["water_year"].values:
        v = pyr["runoff_onset"].sel(water_year=wy).notnull().astype("uint8").rio.write_crs("EPSG:4326")
        valid.append(v.rio.reproject_match(template, resampling=re_.Resampling.nearest).astype(bool)
                     .expand_dims(water_year=[int(wy)]))
    valid = xr.concat(valid, dim="water_year")
    return snow_frac.drop_vars("spatial_ref", errors="ignore"), valid.drop_vars("spatial_ref", errors="ignore")


def zonal_anomalies(config, units_gdf, id_col, unit_dim=None, variables=None,
                    subcells=5, min_weight=0.25, geometry_fixes=None, skip_names=(),
                    anomaly_ds=None, progress=print):
    """Weighted zonal means of the monthly ERA5-Land anomalies per unit
    polygon: weight = polygon coverage fraction x seasonal-snow fraction x
    onset validity (per water year), i.e. the old per-range clip+mask
    semantics as ERA5-cell weights, for ALL units in one pass over the
    store (streamed per variable x 3-month chunk, ~1 GB peak). Cells with
    total weight < ``min_weight`` (in ERA5 cells) are NaN.

    Returns a Dataset dims (``unit_dim`` = ``id_col``, water_year, month)
    with one variable per ERA5 variable plus ``zonal_weight``."""
    import scipy.sparse as sp

    unit_dim = unit_dim or id_col
    if anomaly_ds is None:
        anomaly_ds = open_anomaly(config)
    variables = list(variables or [v for v in anomaly_ds.data_vars if 'month' in anomaly_ds[v].dims])
    progress(f"coverage matrix for {len(units_gdf)} units ...")
    W, ids = unit_coverage_matrix(units_gdf, id_col, anomaly_ds, subcells=subcells,
                                  geometry_fixes=geometry_fixes, skip_names=skip_names)
    progress(f"  {W.shape[0]} units, {W.nnz:,} unit-cell weights")
    snow_frac, valid = pyramid_masks(config, anomaly_ds)
    years = [int(y) for y in anomaly_ds["water_year"].values]
    months = list(anomaly_ds["month"].values)
    # per-year cell weights (seasonal-snow fraction x onset validity); the
    # unit x cell matrix is re-scaled per year on the fly (holding 11 scaled
    # copies would be ~2 GB for the ~4,700 basins)
    snow = snow_frac.values.ravel().astype("float32")
    cellw = {}
    for y in years:
        if y in valid["water_year"].values:
            cellw[y] = snow * valid.sel(water_year=y).values.ravel().astype("float32")
        else:  # a year without onset data yet: seasonal-snow weight only
            cellw[y] = snow

    def weighted(y):
        return sp.csr_matrix(W.multiply(cellw[y][np.newaxis, :]))

    weight_sum = np.stack([np.asarray(weighted(y).sum(axis=1)).ravel() for y in years], axis=1).astype("float32")

    out = {v: np.full((len(ids), len(years), len(months)), np.nan, dtype="float32") for v in variables}
    for v in variables:
        for m0 in range(0, len(months), 3):   # store chunks are 3 months wide
            slab = anomaly_ds[v].isel(month=slice(m0, m0 + 3)).load()
            for j, y in enumerate(years):
                Wy = weighted(y)
                for mi in range(slab.sizes["month"]):
                    x = slab.isel(water_year=j, month=mi).values.ravel()
                    finite = np.isfinite(x)
                    num = Wy @ np.where(finite, x, 0.0)
                    den = Wy @ finite.astype("float64")
                    with np.errstate(invalid="ignore", divide="ignore"):
                        mean = num / den
                    mean[den < min_weight] = np.nan
                    out[v][:, j, m0 + mi] = mean
            del slab
            progress(f"  {v} months {m0}-{m0 + 2} done")
    ds = xr.Dataset({v: ((unit_dim, "water_year", "month"), out[v]) for v in variables},
                    coords={unit_dim: ids, "water_year": years, "month": months})
    ds["zonal_weight"] = ((unit_dim, "water_year"), weight_sum)
    ds["zonal_weight"].attrs["long_name"] = "sum of ERA5-cell weights (coverage x seasonal-snow fraction x onset validity)"
    for v in variables:
        ds[v].attrs.update(anomaly_ds[v].attrs)
        ds[v].attrs["long_name"] = f"{v} anomaly, weighted zonal mean"
    ds.attrs.update({
        "method": ("weighted zonal mean of the monthly ERA5-Land anomaly store; weight = polygon "
                   f"coverage fraction ({subcells}x{subcells} sub-cells) x seasonal-snow fraction "
                   "(pyramid seasonal_snow_pct, level 7) x onset validity (pyramid runoff_onset notnull)"),
        "min_weight_cells": min_weight,
        "anomaly_base_water_years": str(anomaly_ds.attrs.get("anomaly_base_water_years")),
        "dataset_version": config.version, "unit_id_column": id_col,
    })
    return ds


# ---------------------------------------------------------------------------
# on-the-fly Equal Earth derivation (replaces the retired standing stores)

def to_equal_earth(ds):
    """Reproject to EPSG:8857 (Equal Earth) — equal-area pixels for global
    maps/statistics. Derive, don't materialize."""
    import odc.geo.xr  # noqa: F401  -- registers the .odc accessor
    return ds.odc.reproject("EPSG:8857")


def pixelwise_correlations(config, variables=None, level=7, out_path=None,
                           progress=print):
    """Per-pixel Pearson correlation across water years between the
    runoff-onset anomaly and each ERA5-Land variable's monthly anomaly, on
    the Equal-Earth grid of :func:`combined_anomaly_eqearth` and with its
    exact semantics (ERA5 reprojected to EPSG:8857 nearest-neighbour and
    masked per year to onset-valid pixels; every anomaly vs the median over
    all water years; the onset anomaly masked to valid-median pixels; the
    pyramid's yearly ``temporal_resolution`` correlated too and broadcast
    over month) — but streamed one (variable, water year) slab at a time
    (~3 GB peak, no cluster) instead of one dask graph over the whole ~20 GB
    stack, which never converged on a 16 GB box (2026-09-03).

    Writes a zarr with dims ``(month, y, x)`` at ``out_path`` (default
    :func:`correlations_store_path`, local under scratch/) and returns the
    path. The exploratory ``analyses/climate/pixelwise_climate_correlations``
    notebook is its only caller."""
    import odc.geo.xr  # noqa: F401  -- registers the .odc accessor
    import rioxarray  # noqa: F401
    from gsro_analysis.datacube import open_coarse_onset

    stack = open_era5_land(config)
    variables = list(variables or [v for v in VARIABLES if v in stack])
    years = [int(y) for y in stack.water_year.values]
    months = list(stack.month.values)

    # the target grid: what to_equal_earth gives the whole stack (one slab is enough to fix it)
    template = to_equal_earth(stack[variables[0]].isel(water_year=0, month=0)
                              .to_dataset()).compute()
    geobox = template.odc.geobox
    ny, nx = geobox.shape

    # the coarse onset (public pyramid, ~10 km) on that grid, years aligned with the stack
    coarse = open_coarse_onset(config, level)[['runoff_onset', 'runoff_onset_median',
                                                'temporal_resolution']]
    onset = coarse.rio.reproject_match(template).reindex(water_year=years)
    valid = onset['runoff_onset'].notnull()                          # (water_year, y, x)
    onset_anom = (onset['runoff_onset'] - onset['runoff_onset'].median('water_year')
                  ).where(onset['runoff_onset_median'] > 0)
    coords = {'water_year': years, 'y': template.y, 'x': template.x}

    def corr_with_onset(da):
        anom = da - da.median('water_year')
        return xr.corr(onset_anom, anom, dim='water_year').values.astype('float32')

    out = {}
    for v in variables:
        cube = np.empty((len(years), len(months), ny, nx), dtype='float32')
        for j, y in enumerate(years):
            slab = stack[v].sel(water_year=y).load()                 # (month, lat, lon), ~230 MB
            cube[j] = slab.odc.reproject(geobox).values
            del slab
        r = np.full((len(months), ny, nx), np.nan, dtype='float32')
        for mi in range(len(months)):
            da = xr.DataArray(cube[:, mi], dims=('water_year', 'y', 'x'), coords=coords).where(valid)
            r[mi] = corr_with_onset(da)
            del da
        out[v] = (('month', 'y', 'x'), r)
        del cube
        progress(f"{v} done")
    tr = corr_with_onset(onset['temporal_resolution'])
    out['temporal_resolution'] = (('month', 'y', 'x'), np.broadcast_to(tr, (len(months), ny, nx)).copy())

    ds = xr.Dataset(out, coords={'month': months, 'y': template.y, 'x': template.x,
                                 'spatial_ref': template['spatial_ref']})
    ds.attrs.update({
        'method': ('Pearson r across water years of the runoff-onset anomaly (public pyramid level '
                   f'{level}, nearest) vs each ERA5-Land monthly anomaly, both on the Equal Earth grid of '
                   'to_equal_earth; ERA5 masked per year to onset-valid pixels; anomalies vs the median '
                   'over all water years; onset anomaly masked to valid-median pixels'),
        'water_years': years, 'dataset_version': config.version,
    })
    for name in list(ds.variables):
        ds[name].encoding = {}
    ds = ds.chunk({'month': -1, 'y': 512, 'x': 512})
    path = out_path or correlations_store_path(config)
    ds.to_zarr(path, mode='w')
    progress(f"wrote {path}")
    return path


def combined_anomaly_eqearth(config, coarse_onset_ds=None):
    """The on-the-fly replacement for the retired
    ``combined_runoff_onset_and_era5_eqearth_anomaly_ds`` store: raw ERA5
    stack reprojected to Equal Earth, masked to onset-valid pixels, merged
    with the coarse onset dataset, then anomalies vs the all-year median
    (onset anomaly masked to valid-median pixels) — the exact semantics of
    the old combine + anomaly steps, computed lazily.

    ``coarse_onset_ds``: a coarse-resolution onset dataset; defaults to
    ``datacube.open_coarse_onset(config, level=7)`` (~10 km, matches
    ERA5-Land).
    """
    if coarse_onset_ds is None:
        from gsro_analysis.datacube import open_coarse_onset
        coarse_onset_ds = open_coarse_onset(config, level=7)

    era5_proj = to_equal_earth(open_era5_land(config))
    onset_proj = coarse_onset_ds.rio.reproject_match(era5_proj)
    era5_masked = era5_proj.where(onset_proj['runoff_onset'].notnull())
    combined = xr.merge([onset_proj, era5_masked], compat='override')

    static = [v for v in ('runoff_onset_median', 'runoff_onset_mad',
                          'temporal_resolution_median') if v in combined]
    median = combined.drop_vars(static).median(dim='water_year')
    anomaly = combined - median
    if 'runoff_onset_median' in combined:
        anomaly['runoff_onset'] = anomaly['runoff_onset'].where(
            combined['runoff_onset_median'] > 0)
    else:
        warnings.warn("no runoff_onset_median in coarse onset dataset; "
                      "onset anomaly left unmasked")
    return anomaly
