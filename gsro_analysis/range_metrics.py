"""Per-range statistics of the mountain-range cube: the functions behind ``range_statistics.ipynb``, which writes
the wide per-range table (``mountain_range_statistics.csv``), its data dictionary, the per-band hypsometry and the
across-range summary (``mountain_range_summary.json``).

Every function takes the cube (or tables derived from it) as an argument and returns a pandas object indexed by
range name; nothing here reads files. The rules match the per-range metrics table of
``0_aggregate_by_mountain_range.ipynb``, which the wide table extends:

- a (range, elevation, aspect) bin counts only with at least :data:`MIN_PIXELS_PER_BIN` analysed pixels
  (``aggregate.threshold``);
- a range-year counts only if its valid pixels are at least :data:`MIN_YEAR_FRACTION` of the range's
  median-onset pixels (``coverage``); range-years below :data:`LOW_COVERAGE` are flagged, not dropped;
- the per-year topographic fits (lapse rate and aspect amplitude of one water year) additionally keep a bin-year
  only if it holds more than :data:`MIN_BIN_YEAR_FRACTION` of the bin's median-onset pixels, so a bin that is only
  partly observed that year does not bend the year's fit.

Pixel counts are area: the fleet counts pixels of an 80 m UTM grid, so one pixel is :data:`PIXEL_AREA_KM2` at
every latitude (never apply cos(latitude)).
"""

import numpy as np
import pandas as pd
import scipy.stats as sps
import xarray as xr

from gsro_analysis import aggregate

PIXEL_AREA_KM2 = 0.08 ** 2
MIN_PIXELS_PER_BIN = 100
MIN_YEAR_FRACTION = 0.1
LOW_COVERAGE = 0.3
MIN_BIN_YEAR_FRACTION = 0.3
MIN_ASPECTS = 16                 # aspect bins (of 24) a curve needs before a cosine is fitted
MIN_YEARS_FOR_FIT = 3            # paired years for a temperature-sensitivity fit
MIN_YEARS_INTERANNUAL = 8        # years with a value for an interannual s.d. or a co-variation
INCLUSION_SNOW_PCT = 2.0         # a range is analysed only if mapped snow covers at least this share of its polygon
SPRING_MONTHS = ["spring_month_1", "spring_month_2", "spring_month_3"]
BIN_DIMS = ["elevation", "aspect"]
LAT_BANDS = [0, 35, 50, 60, 90]


def base_view(ds, min_pixels=MIN_PIXELS_PER_BIN):
    """The cube with bin means (and stds) of bins below ``min_pixels`` masked: the metrics table's view."""
    return aggregate.threshold(ds, min_pixels)


def _names(ds):
    return pd.Index(ds["mountain_range"].values.astype(str), name="name")


# ---------------------------------------------------------------------------------------------------------------
# area, snow and inclusion

def polygon_area_km2(gmba_gdf, id_col="GMBA_V2_ID"):
    """Polygon area per GMBA id in km² (EPSG:6933, equal area)."""
    return pd.Series(gmba_gdf.to_crs("EPSG:6933").area.values / 1e6, index=gmba_gdf[id_col].values, name="polygon_km2")


def area_totals(hypsometry_df):
    """Per range, the column sums of :func:`hypsometry` (flat pixels included): ``terrain_km2``, ``snow_km2`` (the
    analysis footprint), ``snow_under_forest_km2``, and ``snow_pct_of_terrain``."""
    out = hypsometry_df.groupby("name")[["terrain_km2", "snow_km2", "snow_under_forest_km2"]].sum()
    out["snow_pct_of_terrain"] = 100 * out["snow_km2"] / out["terrain_km2"].where(out["terrain_km2"] > 0)
    return out


def inclusion(snow_km2, polygon_km2, has_metrics, rule_pct=INCLUSION_SNOW_PCT):
    """``snow_pct_of_polygon`` and the inclusion flag: the range has metrics and mapped snow (FCF <= 50 %) covers
    at least ``rule_pct`` % of its GMBA polygon."""
    pct = 100 * snow_km2 / polygon_km2
    return pd.DataFrame({"snow_pct_of_polygon": pct, "included": has_metrics & (pct >= rule_pct)})


def hypsometry(terrain_ds, open_ds, forest_ds):
    """Long table (range, elevation band) of km²: ``terrain_km2`` (land, any forest cover), ``snow_km2`` (analysed,
    FCF <= 50 %) and ``snow_under_forest_km2`` (analysed, FCF > 50 %). Inputs are elevation-only views."""
    def km2(ds, var, label):
        return (ds[var].fillna(0) * PIXEL_AREA_KM2).to_series().rename(label)
    t = km2(terrain_ds, "terrain_n", "terrain_km2")
    s = km2(open_ds, "runoff_onset_median_n", "snow_km2")
    f = km2(forest_ds, "runoff_onset_median_n", "snow_under_forest_km2")
    df = pd.concat([t, s, f], axis=1).fillna(0)
    df = df[(df > 0).any(axis=1)].reset_index().rename(columns={"mountain_range": "name"})
    return df


# ---------------------------------------------------------------------------------------------------------------
# coverage and the range-year anomaly

def coverage(view_ds):
    """Range x water-year share of the range's median-onset pixels with a valid onset that year."""
    num = view_ds["runoff_onset_n"].sum(BIN_DIMS)
    den = view_ds["runoff_onset_median_n"].sum(BIN_DIMS)
    return (num / den.where(den > 0)).to_pandas().set_axis(_names(view_ds))


def range_year_anomaly(view_ds, min_year_fraction=MIN_YEAR_FRACTION):
    """Range x water-year pixel-weighted mean onset anomaly (days), masked where coverage < ``min_year_fraction``:
    the metrics table's anomaly."""
    anom = aggregate.weighted_mean(view_ds, "runoff_onset_anomaly", BIN_DIMS)
    cov = coverage(view_ds)
    return anom.to_pandas().set_axis(_names(view_ds)).where(cov >= min_year_fraction)


def spring_temperature(view_ds, months=SPRING_MONTHS):
    """Range x water-year mean spring 2 m temperature anomaly (°C) from the cube's ERA5-Land zonal means."""
    return view_ds["temperature_2m"].sel(month=months).mean("month").to_pandas().set_axis(_names(view_ds))


def interannual_anomaly(anom_df, min_years=MIN_YEARS_INTERANNUAL):
    """Per range: ``n_years`` with an anomaly, its s.d. and max - min (NaN below ``min_years``), the earliest and the
    latest water year."""
    n = anom_df.notna().sum(axis=1)
    ok = n >= min_years
    out = pd.DataFrame({"n_years": n,
                        "anomaly_sd": anom_df.std(axis=1).where(ok),
                        "anomaly_spread": (anom_df.max(axis=1) - anom_df.min(axis=1)).where(ok)})
    out["earliest_year"] = anom_df[ok].idxmin(axis=1).reindex(out.index)
    out["latest_year"] = anom_df[ok].idxmax(axis=1).reindex(out.index)
    return out


# ---------------------------------------------------------------------------------------------------------------
# aspect: the cosine "sunny-shady difference"

def aspect_curve(onset, n, min_aspects=MIN_ASPECTS):
    """Elevation-relative onset per aspect bin: each band's onset minus the band's median across aspects (bands with
    at least ``min_aspects`` aspect bins), averaged over bands weighted by pixel count. Dims ``(..., elevation,
    aspect)`` -> ``(..., aspect)``."""
    ok = onset.notnull().sum("aspect") >= min_aspects
    rel = (onset - onset.median("aspect")).where(ok)
    w = n.where(rel.notnull())
    return (rel * w).sum("elevation") / w.sum("elevation").where(w.sum("elevation") > 0)


def cosine_fit(curve, aspect_deg, min_aspects=MIN_ASPECTS):
    """Fit ``curve = a + A cos(aspect - phi)``: returns (A [days], azimuth of the EARLIEST onset [deg], r²).
    The sunny-shady difference is 2A."""
    y = np.asarray(curve, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < min_aspects:
        return np.nan, np.nan, np.nan
    theta = np.deg2rad(np.asarray(aspect_deg, dtype=float)[ok])
    X = np.column_stack([np.ones(ok.sum()), np.cos(theta), np.sin(theta)])
    b = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    amp = np.hypot(b[1], b[2])
    latest = np.degrees(np.arctan2(b[2], b[1])) % 360
    resid = y[ok] - X @ b
    r2 = 1 - (resid ** 2).sum() / ((y[ok] - y[ok].mean()) ** 2).sum()
    return amp, (latest + 180) % 360, r2


def yearly_view(view_ds, min_bin_year_fraction=MIN_BIN_YEAR_FRACTION):
    """``runoff_onset`` and ``runoff_onset_n`` with bin-years of no more than ``min_bin_year_fraction`` of the
    bin's median-onset pixels masked (for the per-year topographic fits)."""
    enough = view_ds["runoff_onset_n"] > min_bin_year_fraction * view_ds["runoff_onset_median_n"]
    return view_ds["runoff_onset"].where(enough), view_ds["runoff_onset_n"].where(enough)


def sunny_shady(view_ds):
    """Per range: ``sunny_shady_d`` (2A of the climatological curve, days), ``earliest_aspect_deg``, ``aspect_fit_r2``.
    Range x water year: the same 2A per year, and the year's earliest-onset aspect (two DataFrames)."""
    aspect = view_ds["aspect"].values
    med = aspect_curve(view_ds["runoff_onset_median"], view_ds["runoff_onset_median_n"])
    rows = [cosine_fit(med.sel(mountain_range=r).values, aspect) for r in view_ds["mountain_range"].values]
    clim = pd.DataFrame(rows, columns=["amp", "earliest_aspect_deg", "aspect_fit_r2"], index=_names(view_ds))
    clim.insert(0, "sunny_shady_d", 2 * clim.pop("amp"))
    onset_y, n_y = yearly_view(view_ds)
    yr = aspect_curve(onset_y, n_y)
    years = view_ds["water_year"].values.astype(int)
    fits = [[cosine_fit(yr.sel(mountain_range=r, water_year=y).values, aspect) for y in years]
            for r in view_ds["mountain_range"].values]
    per_year = pd.DataFrame([[2 * f[0] for f in row] for row in fits], index=_names(view_ds), columns=years)
    per_year_az = pd.DataFrame([[f[1] for f in row] for row in fits], index=_names(view_ds), columns=years)
    return clim, per_year, per_year_az


# ---------------------------------------------------------------------------------------------------------------
# lapse rate per water year and its co-variation with the aspect amplitude

def lapse_fit(onset, n, elevation, min_bins=5, min_span_m=500):
    """Count-weighted OLS of the aspect-mean onset on elevation (one range, one field ``(elevation, aspect)``):
    slope in days per 100 m, or NaN with fewer than ``min_bins`` bands or a span below ``min_span_m``."""
    w = n.where(onset.notnull())
    ym = ((onset * w).sum("aspect") / w.sum("aspect")).values
    ww = w.sum("aspect").values
    ok = np.isfinite(ym) & (ww > 0)
    if ok.sum() < min_bins or np.ptp(elevation[ok]) < min_span_m:
        return np.nan
    X = np.column_stack([np.ones(ok.sum()), elevation[ok]])
    W = np.sqrt(ww[ok])
    b = np.linalg.lstsq(X * W[:, None], ym[ok] * W, rcond=None)[0]
    return b[1] * 100


def lapse_rate_weighted_bins(view_ds, min_pixels=MIN_PIXELS_PER_BIN):
    """Per range, the metrics table's lapse rate: sqrt(count)-weighted regression of the bin median onset on elevation
    over every (elevation, aspect) bin with more than ``min_pixels`` pixels (days per 100 m; NaN with fewer than 2
    bins or a span under 100 m)."""
    out = {}
    elevation_axis = view_ds["elevation"].values
    for name in view_ds["mountain_range"].values:
        r = view_ds.sel(mountain_range=name)
        onset = r["runoff_onset_median"].values
        counts = r["runoff_onset_median_n"].values.astype(float)
        elevation = np.broadcast_to(elevation_axis[:, None], onset.shape)
        ok = np.isfinite(onset) & (counts > min_pixels)
        slope = np.nan
        if ok.sum() >= 2 and np.ptp(elevation[ok]) >= 100:
            slope = np.polyfit(elevation[ok], onset[ok], 1, w=np.sqrt(counts[ok]))[0] * 100
        out[str(name)] = slope
    return pd.Series(out, name="lapse_rate_weighted_bins_per_100m").rename_axis("name")


def yearly_topography(view_ds):
    """Range x water year lapse rate (d/100 m) and sunny-shady difference (2A, days), from :func:`yearly_view`."""
    elevation = view_ds["elevation"].values
    onset_y, n_y = yearly_view(view_ds)
    years = view_ds["water_year"].values.astype(int)
    lapse = pd.DataFrame([[lapse_fit(onset_y.sel(mountain_range=r, water_year=y), n_y.sel(mountain_range=r, water_year=y),
                                     elevation) for y in years] for r in view_ds["mountain_range"].values],
                         index=_names(view_ds), columns=years)
    _, amp, _ = sunny_shady(view_ds)
    return lapse, amp


def topography_covariation(lapse_y, amp_y, min_years=MIN_YEARS_INTERANNUAL):
    """Per range: interannual s.d. of the lapse rate and of the sunny-shady difference, and the Pearson r (and p)
    between them across the years where both exist (NaN below ``min_years``)."""
    rows = {}
    for r in lapse_y.index:
        both = lapse_y.loc[r].notna() & amp_y.loc[r].notna()
        n = int(both.sum())
        row = {"topography_n_years": n, "lapse_rate_sd": np.nan, "sunny_shady_sd": np.nan,
               "lapse_aspect_r": np.nan, "lapse_aspect_p": np.nan}
        if n >= min_years:
            x, y = lapse_y.loc[r][both].values, amp_y.loc[r][both].values
            pr = sps.pearsonr(x, y)
            row.update(lapse_rate_sd=x.std(ddof=1), sunny_shady_sd=y.std(ddof=1),
                       lapse_aspect_r=pr.statistic, lapse_aspect_p=pr.pvalue)
        rows[r] = row
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("name")


def covariation_null(lapse_y, amp_y, ranges, n_perm=2000, seed=0, min_years=MIN_YEARS_INTERANNUAL):
    """Shuffle null of the lapse-aspect co-variation: the amplitude years are permuted within each range. Returns the
    observed median r, share positive and number with p < 0.05, and the null's 95th percentiles and p-values."""
    rng = np.random.default_rng(seed)
    pairs = []
    for r in ranges:
        both = lapse_y.loc[r].notna() & amp_y.loc[r].notna()
        if both.sum() >= min_years:
            pairs.append((lapse_y.loc[r][both].values, amp_y.loc[r][both].values))
    n = np.array([len(x) for x, _ in pairs])

    def stats_of(rs):
        t = np.abs(rs) * np.sqrt((n - 2) / np.clip(1 - rs ** 2, 1e-12, None))
        p = 2 * sps.t.sf(t, n - 2)
        return np.median(rs), (rs > 0).mean(), int((p < 0.05).sum())

    obs = stats_of(np.array([np.corrcoef(x, y)[0, 1] for x, y in pairs]))
    null = np.array([stats_of(np.array([np.corrcoef(x, rng.permutation(y))[0, 1] for x, y in pairs]))
                     for _ in range(n_perm)])
    return {"n_ranges": len(pairs), "median_r": obs[0], "frac_positive": obs[1], "n_p_lt_005": obs[2],
            "null_median_r_p95": float(np.percentile(null[:, 0], 95)), "p_median_r": float((null[:, 0] >= obs[0]).mean()),
            "null_frac_positive_p95": float(np.percentile(null[:, 1], 95)), "p_frac_positive": float((null[:, 1] >= obs[1]).mean()),
            "null_n_p_lt_005_p95": float(np.percentile(null[:, 2], 95)), "n_perm": n_perm, "seed": seed}


# ---------------------------------------------------------------------------------------------------------------
# spring-temperature sensitivity

def temperature_sensitivity(anom_df, temp_df, min_years=MIN_YEARS_FOR_FIT, n_perm=10000, seed=0):
    """Per range, the range-year onset anomaly regressed on the spring temperature anomaly over the paired years:
    OLS slope (d/°C) with its standard error and t-based 95 % CI, Pearson r and p; the Theil-Sen slope with its 95 %
    bounds; Spearman rho and p; and a permutation p of |r| (``n_perm`` random permutations of the years)."""
    rng = np.random.default_rng(seed)
    rows = {}
    for r in anom_df.index:
        d = pd.DataFrame({"x": temp_df.loc[r].values, "y": anom_df.loc[r].values}).dropna()
        n = len(d)
        row = dict.fromkeys(["ols_slope", "ols_se", "ols_ci_low", "ols_ci_high", "pearson_r", "ols_p",
                             "theil_sen_slope", "theil_sen_low", "theil_sen_high", "spearman_rho", "spearman_p",
                             "permutation_p"], np.nan)
        row["sensitivity_n"] = n
        if n >= min_years:
            x, y = d["x"].values, d["y"].values
            ols = sps.linregress(x, y)
            tcrit = sps.t.ppf(0.975, n - 2) if n > 2 else np.nan
            ts = sps.mstats.theilslopes(y, x, 0.95)
            sp = sps.spearmanr(x, y)
            perm = np.array([np.corrcoef(rng.permutation(x), y)[0, 1] for _ in range(n_perm)])
            row.update(ols_slope=ols.slope, ols_se=ols.stderr, ols_ci_low=ols.slope - tcrit * ols.stderr,
                       ols_ci_high=ols.slope + tcrit * ols.stderr, pearson_r=ols.rvalue, ols_p=ols.pvalue,
                       theil_sen_slope=ts[0], theil_sen_low=ts[2], theil_sen_high=ts[3],
                       spearman_rho=sp.statistic, spearman_p=sp.pvalue,
                       permutation_p=(np.abs(perm) >= abs(ols.rvalue) - 1e-12).mean())
        rows[r] = row
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("name")


def field_significance(anom_df, temp_df, ranges, n_perm=1000, seed=0, min_years=MIN_YEARS_FOR_FIT):
    """Shared-year permutation null for the count of significantly negative ranges: one permutation of the
    water-year labels of the temperature series, the same for every range (keeps the spatial dependence between
    ranges). Returns the observed count, median slope and share negative, and the null's 95th percentile and max."""
    rng = np.random.default_rng(seed)
    Y = anom_df.loc[ranges]
    X = temp_df.loc[ranges]
    years = list(X.columns)

    def count(Xp):
        neg, slopes = 0, []
        for r in Y.index:
            ok = Y.loc[r].notna().values & Xp.loc[r].notna().values
            if ok.sum() < min_years:
                continue
            fit = sps.linregress(Xp.loc[r].values[ok], Y.loc[r].values[ok])
            slopes.append(fit.slope)
            neg += (fit.pvalue < 0.05) and (fit.slope < 0)
        return neg, np.median(slopes), np.mean(np.array(slopes) < 0)

    obs = count(X)
    null = np.array([count(X.set_axis(list(rng.permutation(years)), axis=1)[years]) for _ in range(n_perm)])
    return {"obs_sig_neg": int(obs[0]), "null_sig_neg_p95": float(np.percentile(null[:, 0], 95)),
            "null_sig_neg_max": int(null[:, 0].max()), "obs_median_slope": float(obs[1]),
            "null_median_slope_min": float(null[:, 1].min()), "obs_frac_neg": float(obs[2]),
            "null_frac_neg_max": float(null[:, 2].max()), "n_perm": n_perm, "seed": seed}


# ---------------------------------------------------------------------------------------------------------------
# Evan & Eisenman (2021): the expected sensitivity from the annual temperature cycle

EE_OMEGA = 2 * np.pi / 365.0
EE_RAM = 0.34
EE_TM = 0.18


def evan_eisenman_fields(clim_ds, tm=EE_TM, ram=EE_RAM):
    """The two expectations of the Evan & Eisenman (2021) idealized model on the T0 / T1 climatology grid (days per °C,
    NaN where the annual cycle never crosses ``tm``). The model: T(t) = T0 - T1 sin(wt - phi); snow accumulates at
    Ra while T < Tm (until precipitation stops) and melts at Rm while T > Tm.

    - ``disappearance``: their eq. 2, d(zeta)/dT0 = -(1/w)(1 + Ra/Rm) / sqrt(T1^2 - (T0 - Tm)^2), the snow-disappearance
      date; the factor (1 + Ra/Rm) is the shorter accumulation season (less snow to melt).
    - ``melt_onset``: the start of melt, the upward crossing of Tm, t = (pi + phi - arcsin((T0 - Tm)/T1)) / w, so
      d(t)/dT0 = -(1/w) / sqrt(T1^2 - (T0 - Tm)^2): eq. 2 without the accumulation factor, independent of Ra/Rm. This is
      the like-for-like expectation for runoff onset. (Both forms verified by simulating the model and differencing.)"""
    graze = clim_ds["T1"] ** 2 - (clim_ds["T0"] - tm) ** 2
    onset = (-(1 / EE_OMEGA) / np.sqrt(graze.where(graze > 0)))
    return onset * (1 + ram), onset


def evan_eisenman_expected(geometries, clim_ds, min_cells=3, tm=EE_TM, ram=EE_RAM):
    """Per range (``geometries``: a GeoSeries in EPSG:4326 indexed by range name), the median over the climatology cells
    the range polygon touches of both Evan & Eisenman (2021) expectations (:func:`evan_eisenman_fields`): ``ee_n_cells``
    (cells where the cycle crosses the threshold), ``ee_defined`` (at least ``min_cells``), ``ee_expected_sensitivity``
    (eq. 2, snow disappearance) and ``ee_expected_onset_sensitivity`` (melt onset), days per °C. The median, not the
    mean, so cells near the eq.'s singularity (T1 ~ |T0 - Tm|) do not dominate."""
    import rasterio.features
    import rasterio.transform
    disappearance, onset = evan_eisenman_fields(clim_ds, tm, ram)
    lat, lon = clim_ds["latitude"].values, clim_ds["longitude"].values
    if lon.max() > 180:
        raise ValueError("expected longitudes in -180..180")
    flip = lat[0] < lat[-1]                                    # rasterize needs north-up rows
    on = onset.values[::-1] if flip else onset.values
    dlon, dlat = abs(lon[1] - lon[0]), abs(lat[1] - lat[0])
    transform = rasterio.transform.from_origin(lon.min() - dlon / 2, lat.max() + dlat / 2, dlon, dlat)
    rows = {}
    for name, geom in geometries.items():
        mask = rasterio.features.rasterize([(geom, 1)], out_shape=on.shape, transform=transform, all_touched=True).astype(bool)
        v = on[mask]
        v = v[np.isfinite(v)]
        ok = len(v) >= min_cells
        med = float(np.median(v)) if ok else np.nan
        rows[name] = {"ee_n_cells": int(len(v)), "ee_defined": ok,
                      "ee_expected_sensitivity": med * (1 + ram), "ee_expected_onset_sensitivity": med}
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("name")


# ---------------------------------------------------------------------------------------------------------------
# the seasonal arch: aspect amplitude per elevation band, by the band's melt date

def band_aspect(view_ds):
    """Long table (range, water year or 'median', elevation band): the band's pixel count, aspect-mean onset and
    cosine amplitude A / earliest azimuth of the band's own aspect profile (elevation-relative), for bands with at
    least :data:`MIN_ASPECTS` aspect bins. Per-year rows use :func:`yearly_view`."""
    aspect, elev = view_ds["aspect"].values, view_ds["elevation"].values
    onset_y, n_y = yearly_view(view_ds)

    def per_band(onset, n):
        rel = (onset - onset.median("aspect")).values
        w = n.where(onset.notnull())
        mean = ((onset * w).sum("aspect") / w.sum("aspect")).values
        enough = (onset.notnull().sum("aspect") >= MIN_ASPECTS).values
        out = []
        for i in np.nonzero(enough)[0]:
            amp, az, _ = cosine_fit(rel[i], aspect)
            out.append((elev[i], mean[i], amp, az))
        return out

    rows = []
    for r in view_ds["mountain_range"].values:
        med = view_ds["runoff_onset_median"].sel(mountain_range=r)
        n_med = view_ds["runoff_onset_median_n"].sel(mountain_range=r)
        px = n_med.where(med.notnull()).sum("aspect").to_series()
        for e, m, a, z in per_band(med, n_med):
            rows.append((str(r), "median", e, px[e], m, a, z))
        for y in view_ds["water_year"].values.astype(int):
            for e, m, a, z in per_band(onset_y.sel(mountain_range=r, water_year=y), n_y.sel(mountain_range=r, water_year=y)):
                rows.append((str(r), str(y), e, px[e], m, a, z))
    return pd.DataFrame(rows, columns=["name", "water_year", "elevation", "median_pixels", "onset_dowy", "amp_d",
                                       "earliest_aspect_deg"])


# ---------------------------------------------------------------------------------------------------------------
# forested vs open pixels

def matched_difference(forest_ds, open_ds, var, min_pixels=MIN_PIXELS_PER_BIN):
    """Per range: the pixel-weighted mean of ``forest - open`` over (elevation, aspect) bins where both have at
    least ``min_pixels`` pixels (weight = the smaller count), and the number of matched bins."""
    n_f, n_o = aggregate.counts_for(forest_ds, var), aggregate.counts_for(open_ds, var)
    ok = (n_f >= min_pixels) & (n_o >= min_pixels)
    w = np.minimum(n_f, n_o).where(ok)
    d = (forest_ds[var] - open_ds[var]).where(ok)
    mean = (d * w).sum(BIN_DIMS) / w.sum(BIN_DIMS).where(w.sum(BIN_DIMS) > 0)
    return mean.to_pandas().set_axis(_names(open_ds)), ok.sum(BIN_DIMS).to_pandas().set_axis(_names(open_ds))


def by_latitude_band(abs_lat, bands=LAT_BANDS):
    """Absolute latitude binned into the reporting bands (labels like ``'35-50'``)."""
    return pd.cut(abs_lat, bands, labels=[f"{a}-{b}" for a, b in zip(bands[:-1], bands[1:])])
