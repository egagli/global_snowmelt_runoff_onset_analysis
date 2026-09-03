"""Per-range statistics on the unified aggregate schema (see
``aggregate.reduce_partials``): elevation lapse rates, spring-temperature
sensitivity, anomaly tables. Used by pipeline/scripts/range_metrics.py (the
one durable per-range table, ``results/<version>/mountain_range_metrics.csv``)
and by the mountain_ranges notebooks; :func:`range_metrics_gdf` joins that
table to the GMBA polygons on read (there is no stats geojson product).

Two lapse-rate definitions exist because two analyses use them:

- :func:`lapse_rate_weighted_bins` — sqrt(count)-weighted regression over
  all (elevation, aspect) bins with >= ``min_count`` pixels (the per-range
  lapse-rate bar chart; the pre-2026 ``calculate_mountain_range_lapse_rates``).
- :func:`lapse_rate_profile` — regression of the aspect-collapsed elevation
  profile (the per-range summary table and the lapse-rate choropleth of the
  polar-triplet world map; the pre-2026 era5_analysis cell).
"""

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats as sps

from gsro_analysis import aggregate

SPRING_MONTHS = ["spring_month_1", "spring_month_2", "spring_month_3"]


METRICS_TABLE = "mountain_range_metrics.csv"
GMBA_KEEP = ["GMBA_V2_ID", "MapName", "Level_04", "Hier_Lvl", "Area", "Perimeter", "geometry"]


def range_metrics_table(version):
    """The per-range metrics table ``range_metrics.py`` writes:
    ``analyses/mountain_ranges/results/<version>/mountain_range_metrics.csv``
    (one row per range with data; provenance columns ``_*``)."""
    from gsro_analysis import paths
    path = paths.resultsdir("mountain_ranges", version) / METRICS_TABLE
    if not path.exists():
        raise FileNotFoundError(f"{path} — run pipeline/scripts/range_metrics.py for {version}")
    return pd.read_csv(path)


def range_metrics_gdf(version, gmba_gdf=None):
    """GMBA polygons joined (inner, on ``GMBA_V2_ID``) with the per-range
    metrics table of ``version`` — the choropleth input of the world maps
    and the per-range comparison notebooks. Replaces the stats geojson
    product (dropped 2026-09-03): the join is a second, not a file."""
    from gsro_analysis import aggregate
    gmba = aggregate.load_gmba() if gmba_gdf is None else gmba_gdf
    keep = [c for c in GMBA_KEEP if c in gmba.columns]
    metrics = range_metrics_table(version)
    drop = [c for c in metrics.columns if c in keep and c != "GMBA_V2_ID"]
    return gmba[keep].merge(metrics.drop(columns=drop), on="GMBA_V2_ID", how="inner")


def range_mean_anomaly(mountains_ds, min_year_fraction=None):
    """Pixel-weighted mean onset anomaly per range and water year
    (``runoff_onset_mean_anomaly`` in the old notebooks). With
    ``min_year_fraction`` years whose valid-pixel count is below that
    fraction of the range's median-pixel count are masked (0.1 = the
    analyses' rule)."""
    dims = [d for d in ("elevation", "aspect", "chili_class") if d in mountains_ds.dims]
    anom = aggregate.weighted_mean(mountains_ds, "runoff_onset_anomaly", dims)
    if min_year_fraction is not None:
        frac = (mountains_ds["runoff_onset_n"].sum(dims)
                / mountains_ds["runoff_onset_median_n"].sum(dims))
        anom = anom.where(frac >= min_year_fraction)
    return anom


def build_anomalies_df(mountains_ds):
    """Per-range annual anomaly table (columns = water years, from the
    dataset — never hardcoded) with continent, centroid latitude, and a
    ``mad`` column (median |anomaly| across years)."""
    anom = range_mean_anomaly(mountains_ds)
    water_years = list(anom.water_year.values)
    df = pd.DataFrame(anom.values, index=anom.mountain_range.values, columns=water_years)
    df["continent"] = mountains_ds["continent"].values
    df["latitude"] = mountains_ds["centroid_latitude"].values
    df["name"] = anom.mountain_range.values
    df["mad"] = df[water_years].abs().median(axis=1)
    return df


def lapse_rate_weighted_bins(mountains_ds, min_count=100, min_relief=100):
    """Weighted (sqrt-count) linear onset-vs-elevation regression per
    range over the elevation x aspect bins (CHILI collapsed): lapse rate in
    days per 100 m + weighted R^2. Requires >= 2 bins with count > min_count
    spanning >= min_relief m."""
    ds = aggregate.collapse(mountains_ds) if "chili_class" in mountains_ds.dims else mountains_ds
    rows = []
    for rng in ds.mountain_range.values:
        d = ds.sel(mountain_range=rng)
        onset = d["runoff_onset_median"].values
        counts = d["runoff_onset_median_n"].values.astype(float)
        elev = np.broadcast_to(d["elevation"].values[:, None], onset.shape)
        ok = np.isfinite(onset) & (counts > min_count)
        lapse, r2 = np.nan, np.nan
        if ok.sum() >= 2 and (elev[ok].max() - elev[ok].min()) >= min_relief:
            w = np.sqrt(counts[ok])
            coeffs = np.polyfit(elev[ok], onset[ok], 1, w=w)
            pred = np.polyval(coeffs, elev[ok])
            mean_w = np.average(onset[ok], weights=w)
            ss_tot = np.sum(w * (onset[ok] - mean_w) ** 2)
            ss_res = np.sum(w * (onset[ok] - pred) ** 2)
            r2 = 1 - ss_res / ss_tot
            lapse = coeffs[0] * 100
        rows.append({"name": rng, "lapse_rate": lapse, "r_squared": r2,
                     "CONTINENT": str(d["continent"].values),
                     "latitude": float(d["centroid_latitude"]), "longitude": float(d["centroid_longitude"])})
    return pd.DataFrame(rows)


def elevation_profile(mountains_ds, var="runoff_onset_median"):
    """Aspect- (and CHILI-) collapsed per-elevation mean per range."""
    dims = [d for d in ("aspect", "chili_class") if d in mountains_ds.dims]
    return aggregate.weighted_mean(mountains_ds, var, dims)


def lapse_rate_profile(mountains_ds, min_bins=3):
    """Per-range linear regression of the elevation profile (days/100 m),
    Pearson r and the number of elevation bins used."""
    prof = elevation_profile(mountains_ds)
    rows = []
    for rng in prof.mountain_range.values:
        y = prof.sel(mountain_range=rng).values
        ok = np.isfinite(y)
        if ok.sum() < min_bins:
            rows.append({"name": rng, "lapse_rate_per_100m": np.nan, "correlation": np.nan, "n": int(ok.sum())})
            continue
        res = sps.linregress(prof.elevation.values[ok], y[ok])
        rows.append({"name": rng, "lapse_rate_per_100m": res.slope * 100, "correlation": res.rvalue, "n": int(ok.sum())})
    return pd.DataFrame(rows).set_index("name")


def seasonal_means(climate_ds):
    """Add ``winter/spring/summer_months_mean`` entries on the month axis."""
    parts = [climate_ds]
    for season in ("winter", "spring", "summer"):
        months = [f"{season}_month_{i}" for i in (1, 2, 3)]
        parts.append(climate_ds.sel(month=months).mean("month")
                     .expand_dims(month=[f"{season}_months_mean"]))
    return xr.concat(parts, dim="month")


def spring_temperature_sensitivity(mountains_ds, var="temperature_2m",
                                   months=SPRING_MONTHS, min_years=3,
                                   min_year_fraction=0.1):
    """Per-range regression of the annual mean onset anomaly on the mean
    spring anomaly of ``var`` (default 2 m temperature): OLS slope
    (days per unit), Pearson r, p-value, n, plus the Theil-Sen slope and
    its 95 % bounds. Years with < ``min_year_fraction`` of the range's
    pixels valid are excluded (the analyses' rule)."""
    anom = range_mean_anomaly(mountains_ds, min_year_fraction=min_year_fraction)
    clim = mountains_ds[var].sel(month=list(months)).mean("month")
    rows = []
    for rng in mountains_ds.mountain_range.values:
        df = pd.DataFrame({"x": clim.sel(mountain_range=rng).values,
                           "y": anom.sel(mountain_range=rng).values}).dropna()
        row = {"name": rng, "anomaly_n": len(df), "anomaly_slope": np.nan, "anomaly_corr": np.nan,
               "anomaly_pval": np.nan, "theil_sen_slope": np.nan, "theil_sen_low": np.nan,
               "theil_sen_high": np.nan}
        if len(df) >= min_years:
            res = sps.linregress(df["x"], df["y"])
            ts = sps.mstats.theilslopes(df["y"], df["x"], 0.95)
            row.update(anomaly_slope=res.slope, anomaly_corr=res.rvalue, anomaly_pval=res.pvalue,
                       theil_sen_slope=ts[0], theil_sen_low=ts[2], theil_sen_high=ts[3])
        rows.append(row)
    return pd.DataFrame(rows).set_index("name")


def climate_regressions(mountains_ds, climate_vars, min_years=3, min_year_fraction=0.1):
    """Per (range, ERA5 variable, month): n, slope, corr, pval of the annual
    onset anomaly on the variable's anomaly — the heatmap input of the
    temperature_sensitivity notebook."""
    anom = range_mean_anomaly(mountains_ds, min_year_fraction=min_year_fraction)
    params = ["n", "slope", "corr", "pval"]
    ranges = mountains_ds.mountain_range.values
    months = mountains_ds.month.values
    data = {}
    for var in climate_vars:
        arr = np.full((len(ranges), len(months), len(params)), np.nan)
        for i, rng in enumerate(ranges):
            y = anom.sel(mountain_range=rng).values
            for j, month in enumerate(months):
                x = mountains_ds[var].sel(mountain_range=rng, month=month).values
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() >= min_years:
                    res = sps.linregress(x[ok], y[ok])
                    arr[i, j] = [ok.sum(), res.slope, res.rvalue, res.pvalue]
        data[var] = (("mountain_range", "month", "param"), arr)
    return xr.Dataset(data, coords={"mountain_range": ranges, "month": months, "param": params,
                                    "centroid_latitude": mountains_ds["centroid_latitude"],
                                    "centroid_longitude": mountains_ds["centroid_longitude"],
                                    "continent": mountains_ds["continent"]})


# ---------------------------------------------------------------------------
# notebook-facing preparation of the two cubes (one call replaces the 3-4
# threshold/derivation cells every pre-2026 notebook started with)

def prepare_mountain_ranges(mountains_full, min_pixels=100, annual_fraction=0.3,
                            tropical_andes_rule=True):
    """The mountain-range cube as the analyses use it: CHILI collapsed
    (elevation x aspect), bins with <= ``min_pixels`` masked, a bin-year
    kept only if it has > ``annual_fraction`` of the bin's median-pixel
    count, ``runoff_onset_elev_relative`` (deviation from the per-elevation
    median across aspects, masked where any aspect is missing) and
    ``runoff_onset_mean_anomaly`` (pixel-weighted range mean) added, ranges
    without data dropped. ``tropical_andes_rule``: mask tropical-Andes bins
    with median onset >= 250 DOWY below 5000 m (South American ranges with
    centroid latitude > -20), where late-season values are artefacts."""
    ds = aggregate.collapse(mountains_full) if "chili_class" in mountains_full.dims else mountains_full.copy()
    ds = aggregate.threshold(ds, min_pixels + 1)
    annual_thresh = annual_fraction * ds["runoff_onset_median_n"]
    for v in ("runoff_onset", "runoff_onset_anomaly", "runoff_onset_std", "runoff_onset_anomaly_std"):
        ds[v] = ds[v].where(ds["runoff_onset_n"] > annual_thresh)
    if tropical_andes_rule:
        trop = (ds["continent"] == "South America") & (ds["centroid_latitude"] > -20)
        ok = (ds["runoff_onset_median"] < 250) | (ds["elevation"] > 5000) | ~trop
        for v in ("runoff_onset_median", "runoff_onset_median_std", "runoff_onset_mad", "runoff_onset_mad_std",
                  "runoff_onset", "runoff_onset_anomaly", "runoff_onset_std", "runoff_onset_anomaly_std"):
            ds[v] = ds[v].where(ok)
    rel = aggregate.elevation_relative(ds["runoff_onset_median"])
    ds["runoff_onset_elev_relative"] = rel.where(~rel.isnull().any("aspect"))
    ds["runoff_onset_elev_relative"].attrs = {"units": "days",
                                              "long_name": "median onset minus the per-elevation median across aspects"}
    ds["runoff_onset_mean_anomaly"] = range_mean_anomaly(ds)
    ds["runoff_onset_mean_anomaly"].attrs = {"units": "days", "long_name": "pixel-weighted mean onset anomaly per range"}
    has_data = ds["runoff_onset_median"].notnull().any(("elevation", "aspect"))
    return ds.sel(mountain_range=has_data)


def basin_summary(basins_ds, basins_gdf, populations_gdf=None, pixel_area_km2=(80 / 1000) ** 2,
                  min_area_pct=5, min_year_area_pct=1):
    """One GeoDataFrame row per river basin: geometry (dissolved per
    PFAF_ID), population, pixel-weighted means of median onset / MAD /
    yearly onset and anomaly, pixel counts, mapped area and its share of
    the basin (``area_pct``). Means are masked where the mapped share is
    below ``min_area_pct`` (yearly ones below ``min_year_area_pct``) — the
    rule the basin maps have always used."""
    import geopandas as gpd  # noqa: F401
    dims = [d for d in ("elevation", "aspect", "chili_class") if d in basins_ds.dims]
    ids = basins_ds["river_basin"].values
    gdf = basins_gdf[basins_gdf["PFAF_ID"].isin(ids)].copy()
    if populations_gdf is not None and "total_population" in populations_gdf:
        gdf = gdf.merge(populations_gdf[["HYBAS_ID", "total_population"]], on="HYBAS_ID", how="left")
        gdf = gdf.rename(columns={"total_population": "POPULATION"})
    agg = {"SUB_AREA": "sum"}
    if "POPULATION" in gdf:
        agg["POPULATION"] = "sum"
    for c in ("HYBAS_ID", "MAIN_BAS", "ORDER", "ENDO", "COAST"):
        if c in gdf:
            agg[c] = "first"
    gdf = gdf.dissolve(by="PFAF_ID", aggfunc=agg).reset_index().sort_values("PFAF_ID").set_index("PFAF_ID")
    ds = basins_ds.sel(river_basin=gdf.index.values)

    def mean(var):
        return aggregate.weighted_mean(ds, var, dims)

    n_med = ds["runoff_onset_median_n"].sum(dims)
    gdf["pixel_count"] = n_med.values
    gdf["area_pct"] = 100 * n_med.values * pixel_area_km2 / gdf["SUB_AREA"]
    ok = gdf["area_pct"] > min_area_pct
    gdf["runoff_onset_median"] = mean("runoff_onset_median").values
    gdf["runoff_onset_mad"] = mean("runoff_onset_mad").values
    gdf.loc[~ok, ["runoff_onset_median", "runoff_onset_mad"]] = np.nan
    on, an = mean("runoff_onset"), mean("runoff_onset_anomaly")
    n_y = ds["runoff_onset_n"].sum(dims)
    for y in ds.water_year.values:
        pct = 100 * n_y.sel(water_year=y).values * pixel_area_km2 / gdf["SUB_AREA"]
        oky = pct > min_year_area_pct
        gdf[f"runoff_onset_WY{y}"] = np.where(oky, on.sel(water_year=y).values, np.nan)
        gdf[f"runoff_onset_anomaly_WY{y}"] = np.where(oky, an.sel(water_year=y).values, np.nan)
        gdf[f"pixel_count_WY{y}"] = n_y.sel(water_year=y).values
    return gdf.reset_index()
