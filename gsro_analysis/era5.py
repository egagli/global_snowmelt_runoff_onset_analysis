"""ERA5-Land acquisition and derived products — the minimal set.

Exactly two standing stores per dataset version (on Azure under
``settings.ERA5_DATA_PREFIX/<version>/``), plus the tiny WeatherBench2
climatology cache the temperature_sensitivity notebook keeps locally:

1. ``era5_water_year_<yr>.zarr`` — the canonical acquisition: native 0.1°
   lat/lon, monthly, 9 hemisphere-aware months per water year
   (:data:`MONTH_LABELS`), 8 variables (:data:`VARIABLES`).
2. ``era5_land_anomaly_ds.zarr`` — monthly anomalies vs the median over a
   recorded base period; feeds the per-range/per-basin temperature
   sensitivity via :func:`zonal_anomalies`.

Everything else that used to be materialized (the Equal-Earth family:
``combined_runoff_onset_and_era5_eqearth{,_anomaly}_ds``, the local 10-yr
eqearth stacks, the mountain-masked anomaly variant, the 6-month local
per-WY mirrors) is retired — derive on the fly with
:func:`combined_anomaly_eqearth` / :func:`to_equal_earth` instead
(2026-08-24 decision). Layouts are versioned; code assumes >= v10.

Acquisition needs Earth Engine (``settings.initialize_earthengine()``,
service-account key — works headless) and the ``xee < 0.1`` pin that
easysnowdata's ``get_era5`` requires.
"""

import warnings

import easysnowdata
import numpy as np
import xarray as xr

from gsro_analysis import settings

# 9 hemisphere-aware months per water year: NH Dec(wy-1)-Aug(wy),
# SH Jun(wy)-Feb(wy+1), both labeled winter -> spring -> summer.
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


# ---------------------------------------------------------------------------
# paths (versioned under ERA5_DATA_PREFIX/<version>/)

def water_year_store_path(config, water_year):
    return (f"{settings.ERA5_DATA_PREFIX}/{config.version}/"
            f"era5_water_year_{water_year}.zarr")


def anomaly_store_path(config):
    return (f"{settings.ERA5_DATA_PREFIX}/{config.version}/"
            f"era5_land_anomaly_ds.zarr")


def correlations_store_path(config):
    """LOCAL output of the exploratory pixel-wise correlation notebook
    (analyses/climate/pixelwise_climate_correlations.ipynb): a zarr under the
    gitignored scratch/ folder, keyed by dataset version. Its own product,
    read by nothing else, so it never goes to Azure (2026-09-03)."""
    from gsro_analysis import paths
    paths.SCRATCH.mkdir(parents=True, exist_ok=True)
    return paths.SCRATCH / f"{config.version}_runoff_onset_and_era5_eqearth_anomaly_correlations.zarr"


# ---------------------------------------------------------------------------
# acquisition (Earth Engine)

def fetch_water_year(config, water_year, variables=VARIABLES):
    """One water year of monthly ERA5-Land on the native 0.1° grid,
    NH + SH merged with hemisphere-aware month windows. Earth Engine must
    already be initialized (settings.initialize_earthengine())."""
    hemis = [((-180, 0, 180, 90), f"{water_year - 1}-12-01", f"{water_year}-08-31"),
             ((-180, -90, 180, 0), f"{water_year}-06-01", f"{water_year + 1}-02-28")]
    parts = []
    for bbox, start, end in hemis:
        ds = easysnowdata.hydroclimatology.get_era5(
            version="ERA5_LAND", bbox_input=bbox, cadence="MONTHLY",
            start_date=start, end_date=end, initialize_ee=False,
            variables=list(variables))
        ds = (ds.expand_dims({'water_year': [water_year]})
              .rename({'time': 'month'})
              .assign_coords(month=MONTH_LABELS))
        parts.append(ds)
    return xr.merge(parts)


def _marker_path(store_path):
    """Sidecar done-marker for an ERA5 zarr: <prefix>/_complete/<name>.json —
    the same verify-then-mark ledger the ancillary tiles use."""
    prefix, name = store_path.rsplit('/', 1)
    return f"{prefix}/_complete/{name}.json"


def _sample_finite_fraction(ds, variables):
    """Finite fraction of one spring slab of the first variable, over the
    northern half (ERA5-Land is land-only: ~0.3-0.4 when data is present,
    exactly 0 for a metadata-only / interrupted store)."""
    import numpy as np
    da = ds[variables[0]]
    if 'water_year' in da.dims:
        da = da.isel(water_year=0)
    slab = da.sel(month='spring_month_1').isel(latitude=slice(0, 900)).values
    return float(np.isfinite(slab).mean())


def verify_and_mark(config, path, variables=VARIABLES, min_finite=0.1):
    """Return True iff the store holds every variable AND a data sample reads
    back finite; write the done-marker when it does (self-heals stores that
    predate the ledger). Reads through a FRESH fs (see settings.fresh_blob_fs).
    A 0-variable or NaN-only store — what an interrupted write leaves behind
    (2026-08-25, twice) — is NOT complete."""
    import json
    fs = settings.fresh_blob_fs(config)
    marker = _marker_path(path)
    if fs.exists(marker):
        return True
    if not settings.zarr_store_exists(fs, path):
        return False
    try:
        ds = xr.open_zarr(fs.get_mapper(path), chunks=None)
        if not set(variables) <= set(ds.data_vars):
            return False
        frac = _sample_finite_fraction(ds, variables)
    except Exception as e:  # noqa: BLE001 - reported, then treated as not verified
        print(f"verify_and_mark({path}): read failed: {type(e).__name__}: {e}")
        return False
    if frac < min_finite:
        return False
    fs.pipe_file(marker, json.dumps({'store': path, 'variables': sorted(variables),
                                     'sample_finite_fraction': round(frac, 3)}).encode())
    return True


def build_water_year_stores(config, water_years=None, variables=VARIABLES,
                            overwrite=False):
    """Write the per-water-year stores (skip-if-complete; a partial store
    from an interrupted run is rebuilt). The ONE acquisition code path —
    the old notebook had two divergent variants (9-month Azure vs 6-month
    local); this is the 9-month one."""
    settings.initialize_earthengine()
    fs = config.azure_blob_fs
    for wy in (water_years if water_years is not None else config.water_years):
        wy = int(wy)
        path = water_year_store_path(config, wy)
        if not overwrite and verify_and_mark(config, path, variables):
            print(f"skip WY{wy} (complete: {path})")
            continue
        # Fetch one variable at a time (bounds xee's transient peak — the
        # whole-year request was OOM-killed on a 5 GB-free dev box), then
        # write the year ONCE: appending new variables with mode='a' hits
        # zarr v3's ContainsArrayError on the coordinate arrays (ERA5
        # Acquire run 32882895446). Holding the merged year needs a few GB —
        # this is a 16 GB-runner job (era5_acquire.yml), not a laptop one.
        # consolidated=False for the same reason as
        # datacube.save_ancillary_tile (stale-dircache consolidation).
        if settings.zarr_store_exists(fs, path):  # partial, interrupted run
            fs.rm(path, recursive=True)
            fs.invalidate_cache()
        parts = []
        for var in variables:
            print(f"fetching WY{wy} {var} ...")
            parts.append(fetch_water_year(config, wy, [var]))
        print(f"writing WY{wy} ...")
        xr.merge(parts).to_zarr(fs.get_mapper(path), mode='w',
                                consolidated=False)
        del parts
        # verification through a FRESH fs: the shared one can read back an
        # empty group after the EE-heavy fetch (run 32885451802 failed here
        # on a store that was in fact complete)
        if not (verify_and_mark(config, path, variables)
                or settings.verify_in_subprocess(config, 'era5', 'verify_and_mark',
                                                 path, list(variables))):
            raise RuntimeError(f"WY{wy}: post-write verification failed "
                               f"({path})")
        print(f"wrote {path}")


# ---------------------------------------------------------------------------
# opening + derived products

def open_era5_stack(config, chunks='auto'):
    """All water years of the version's ERA5-Land stores, concatenated on
    water_year (native 0.1° grid)."""
    fs = config.azure_blob_fs
    prefix = water_year_store_path(config, 0).rsplit('/', 1)[0]
    files = sorted(fs.glob(f"{prefix}/era5_water_year_*.zarr"))
    if not files:
        raise FileNotFoundError(f"no era5_water_year_*.zarr under {prefix}")
    datasets = [xr.open_zarr(fs.get_mapper(f), chunks=chunks) for f in files]
    return xr.concat(datasets, dim='water_year')


def build_anomaly_store(config, base_water_years=None):
    """The native-grid monthly anomaly store: stack minus the median over
    ``base_water_years`` (default: every year in the stack); the base period
    is recorded in attrs. Memory-bounded on purpose: the median is computed
    one (variable, month) slab at a time (11 years x 1800 x 3600 -> a few
    hundred MB) and held as a small in-memory dataset; the anomaly itself is
    a lazy broadcast subtraction streamed to zarr. dask's own median needs
    every year of a chunk at once and OOM-killed a 5 GB box (2026-08-25).
    Still a ~4 GB job -> run it on a runner (ERA5 Acquire, build_anomaly).
    """
    import numpy as np

    stack = open_era5_stack(config)
    years = [int(y) for y in (base_water_years if base_water_years is not None
                              else stack.water_year.values)]
    base = stack.sel(water_year=years)
    medians = {}
    for var in stack.data_vars:
        slabs = []
        for month in stack.month.values:
            slab = base[var].sel(month=month).load()          # 11 x 1800 x 3600
            slabs.append(slab.median(dim='water_year'))
            del slab
        medians[var] = xr.concat(slabs, dim='month').assign_coords(
            month=stack.month.values)
        print(f"median done: {var}")
    median_ds = xr.Dataset(medians)
    anomaly = stack - median_ds                              # lazy broadcast
    anomaly.attrs['anomaly_base_water_years'] = years
    anomaly.attrs['dataset_version'] = config.version
    for name in list(anomaly.variables):
        anomaly[name].encoding = {}
    anomaly = anomaly.chunk({'water_year': 1, 'month': 3,
                             'latitude': 450, 'longitude': 900})
    fs = settings.fresh_blob_fs(config)
    path = anomaly_store_path(config)
    if settings.zarr_store_exists(fs, path):
        fs.rm(path, recursive=True)
        fs.invalidate_cache()
    anomaly.to_zarr(fs.get_mapper(path), mode='w', consolidated=False)
    variables = list(stack.data_vars)
    if not (verify_and_mark(config, path, variables)
            or settings.verify_in_subprocess(config, 'era5', 'verify_and_mark',
                                             path, variables)):
        raise RuntimeError(f"anomaly store verification failed ({path})")
    print(f"wrote {path}")
    return path


def open_anomaly(config, chunks='auto'):
    """The monthly anomaly store feeding the per-range temperature
    sensitivity; refuses a store that has not passed verify_and_mark."""
    path = anomaly_store_path(config)
    if not verify_and_mark(config, path):
        raise FileNotFoundError(
            f"{path} is missing or unverified (no _complete marker / no "
            "finite data). Run the ERA5 Acquire workflow with build_anomaly, "
            "or era5.build_anomaly_store(config) on a >= 8 GB machine.")
    return xr.open_zarr(config.azure_blob_fs.get_mapper(path),
                        decode_coords='all', chunks=chunks)


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
    variables = list(variables or anomaly_ds.data_vars)
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

    era5_proj = to_equal_earth(open_era5_stack(config))
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
