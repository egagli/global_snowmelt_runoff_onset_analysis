"""The aggregation engine: per-tile pixel tables -> per-tile PARTIAL SUMS
(the fleet "map", :func:`tile_partials`) -> merged per-unit-type cubes (the
"reduce", :func:`reduce_partials`, run by the notebooks
``analyses/<unit>/0_aggregate_by_<unit>.ipynb``), plus the read-side helpers
every analysis notebook uses on the resulting files (:func:`weighted_mean`,
:func:`collapse`, :func:`threshold`, :func:`elevation_relative`).

Design (August 2026 redesign, replacing the parquet-scan stage 2):

- Every statistic the analyses read is *reducible*: bin means and standard
  deviations come from (n, sum, sum of squares), correlations from
  (n, sum x, sum y, sum xy, sum x^2, sum y^2). No analysis ever read the old
  bin-level medians, so they are gone; a fleet job therefore emits ~0.1-1 MB
  of partial sums per tile instead of a 30-50 MB pixel table, and stage 2
  is a laptop-scale reduce (minutes, ~1 GB) rather than a 16 GB dask job.
- ONE schema for every unit type (:data:`UNIT_TYPES`): dims
  ``(unit, <geometric bins>, chili_class, water_year)`` with plain variables
  ``<quantity>`` (bin mean), ``<quantity>_std`` and ``<quantity>_n`` — no
  ``statistic`` axis, no precomputed marginals (``basin_*``,
  ``elev_relative``): marginals are one call to :func:`weighted_mean`.
  The MAP keys every unit type by the finest axes it will ever need
  (:data:`UNIT_TYPES`): mountain ranges by elevation x aspect, river basins
  (HydroBASINS level 6) by elevation, continents by latitude x elevation x
  aspect; all three carry the CHILI insolation class. The REDUCE writes
  output GROUPS (:data:`GROUPS`): the default ``river_basins`` cube is
  level 5 (``PFAF_ID // 10`` — Pfafstetter codes nest by digit prefix, so
  the sums of the level-6 rows ARE the level-5 sums) and the default
  ``continents`` cube is latitude x elevation (aspect summed out);
  ``river_basins_l6`` and ``continents_aspect`` are the finer cubes on demand.
  Decided 2026-09-02, before the first campaign: a coarser axis is a free
  ``groupby`` at reduce time, a finer one is a fleet re-map.
- Pixel filters are named predicate sets (:data:`FILTERS`) applied in the
  map; every partials row and every output records its ``filter_tag``.
  Two tags: ``fcf_lte_50`` (the analyses' rule) and ``full_dataset`` (read
  once, for the share of pixels above 5000 m). The former ``no_trees`` tag
  was read by nothing and left the map on 2026-09-03.
- ``water_years`` is always a parameter (``config.water_years``) — the
  dataset grows a year at a time, never hardcode it.
- No cluster service, no dask: the map runs inside the fleet worker, the
  reduce is pandas + numpy in a notebook. Reading the partials, summing them
  over tiles, naming the units and writing the cubes are visible cells of the
  aggregation notebooks; this module holds only the math.

Semantics that differ from the pre-2026 stage 2 (deliberate, checked
2026-08-26 on the five v10 dry-run tiles): bins are half-open ``[left,
right)`` (numpy convention; the old ``pd.cut`` was ``(left, right]``, so
0 m pixels were dropped and integer elevations exactly on an edge sat one
bin lower); pixels without a CHILI value are class ``none`` (they were
counted as ``neutral``); aspect 360 deg wraps to 0 (north); the FCF
correlation uses only pixels with FCF coverage (the old one mixed in the
-9999 fill); river basins bin elevation only, so flat pixels (undefined
aspect, ~0.07 %) are included where the old aspect binning dropped them.
Continents carry an aspect key since 2026-09-03: flat pixels keep an
undefined (NaN) aspect in the partials, so the default latitude x elevation
cube still counts them and only the on-demand ``continents_aspect`` cube
drops them.
"""

import operator

import numpy as np
import pandas as pd
import xarray as xr

from gsro_analysis import settings

# ---------------------------------------------------------------------------
# named pixel-filter sets (AND-ed predicates on the pixel table).
# 'fcf_lte_50' is the rule every analysis uses; 'full_dataset' is read once.
# The fcf >= 0 bound also excludes the -9999 nodata fill (pixels with no
# forest-cover coverage).
FILTERS = {
    "full_dataset": [
        ("snow_classification", "!=", 4),
        ("esa_worldcover", "not in", [50, 80]),
    ],
    "fcf_lte_50": [
        ("snow_classification", "!=", 4),
        ("esa_worldcover", "not in", [50, 80]),
        ("forest_cover_fraction", ">=", 0),
        ("forest_cover_fraction", "<=", 50),
    ],
}

_OPS = {
    "==": operator.eq, "!=": operator.ne,
    ">": operator.gt, ">=": operator.ge, "<": operator.lt, "<=": operator.le,
    "in": lambda s, v: s.isin(v), "not in": lambda s, v: ~s.isin(v),
}


def apply_filter(df, filter_tag):
    """Rows of a pixel table passing the named filter (``None`` = all rows)."""
    if not filter_tag:
        return df
    mask = pd.Series(True, index=df.index)
    for col, op, val in FILTERS[filter_tag]:
        mask &= _OPS[op](df[col], val)
    return df[mask]


# integer continent encoding in the pixel tables: alphabetical np.unique
# order of the USGS continents shapefile (see datacube.add_mountain_range_
# and_basin_and_continent). The ancillary build stamps this mapping into
# tile attrs; assert against it rather than trusting this copy blindly.
CONTINENTS_ENUM = {
    0: "Africa", 1: "Antarctica", 2: "Asia", 3: "Australia",
    4: "Europe", 5: "North America", 6: "Oceania", 7: "South America",
}
# Australia is folded into Oceania in every output (six continents)
CONTINENT_MERGE = {"Australia": "Oceania"}
CONTINENT_NAMES = ["Africa", "Asia", "Europe", "North America", "Oceania",
                   "South America"]

# CHILI insolation classes (Theobald et al. 2015 empirical thresholds)
CHILI_COOL_MAX = 0.448
CHILI_WARM_MIN = 0.767
CHILI_CLASSES = ["cool", "neutral", "warm", "none"]   # 'none' = no CHILI value

# ---------------------------------------------------------------------------
# bins (edges); coordinates are bin centers
DEM_EDGES = np.arange(0, 9000 + 100, 100)          # 100 m, 0-9000 m
ASPECT_EDGES = np.arange(0, 360 + 15, 15)          # 15 deg
LAT_EDGES = np.arange(-90, 90 + 1, 1)              # 1 deg, whole globe
BIN_EDGES = {"elevation": DEM_EDGES, "aspect": ASPECT_EDGES, "latitude": LAT_EDGES}
BIN_SOURCE = {"elevation": "dem", "aspect": "aspect", "latitude": "original_lat"}
BIN_UNITS = {"elevation": "meters", "aspect": "degrees", "latitude": "degrees"}


def bin_centers(name):
    e = BIN_EDGES[name]
    return e[:-1] + np.diff(e) / 2


# The three unit types the MAP keys partial rows by — the finest axes each
# will ever need. ``bins`` are the geometric axes; a bin in ``optional_bins``
# may be undefined for a pixel (flat pixels have no aspect) without dropping
# the pixel: its row then carries NaN there and only groups that use that
# axis leave it out. ``id_col`` is the ancillary layer: GMBA_V2_ID (standard
# 300 inventory), PFAF_ID at HydroBASINS level :data:`BASIN_LEVEL_STORED`,
# the USGS continent code.
UNIT_TYPES = {
    "mountain_ranges": {"id_col": "GMBA_V2_ID", "unit_dim": "mountain_range",
                        "bins": ("elevation", "aspect"), "optional_bins": ()},
    "river_basins":    {"id_col": "PFAF_ID", "unit_dim": "river_basin",
                        "bins": ("elevation",), "optional_bins": ()},
    "continents":      {"id_col": "continent", "unit_dim": "continent",
                        "bins": ("latitude", "elevation", "aspect"), "optional_bins": ("aspect",)},
}
BASIN_LEVEL_STORED = settings.BASIN_ATLAS_STORED_LEVEL   # 6

# The cubes the REDUCE writes (``reduce_partials(..., group, ...)``): a unit
# type, the bins kept (the others are summed out — exact, the sums are
# additive) and, for basins, the HydroBASINS level (a digit-prefix of the
# stored id). The aggregation notebooks name the file ``all_<group>_<filter_tag>.nc``.
GROUPS = {
    "mountain_ranges":   {"unit_type": "mountain_ranges", "bins": ("elevation", "aspect")},
    "river_basins":      {"unit_type": "river_basins", "bins": ("elevation",), "basin_level": 5},
    "river_basins_l6":   {"unit_type": "river_basins", "bins": ("elevation",), "basin_level": 6},
    "continents":        {"unit_type": "continents", "bins": ("latitude", "elevation")},
    "continents_aspect": {"unit_type": "continents", "bins": ("latitude", "elevation", "aspect")},
}
DEFAULT_GROUPS = ("mountain_ranges", "river_basins", "continents")


def basin_level_ids(pfaf_ids, level, stored=BASIN_LEVEL_STORED):
    """HydroBASINS level-``level`` PFAF_IDs of stored level-``stored`` ids
    (the first ``level`` digits). Every level-k Pfafstetter code has exactly
    k digits, which is asserted."""
    ids = np.asarray(pfaf_ids, dtype="int64")
    if not (1 <= level <= stored):
        raise ValueError(f"basin level {level} not in 1..{stored} (the stored level)")
    lo, hi = 10 ** (stored - 1), 10 ** stored
    if not ((ids >= lo) & (ids < hi)).all():
        raise ValueError(f"PFAF_IDs are not all {stored}-digit level-{stored} codes: "
                         f"min {ids.min()}, max {ids.max()} — partials from an older schema?")
    return ids // 10 ** (stored - level)


NODATA = -9999


# ---------------------------------------------------------------------------
# MAP: one tile's pixel table -> partial sums

def chili_class_index(chili):
    """0 cool / 1 neutral / 2 warm / 3 none, per :data:`CHILI_CLASSES`."""
    idx = np.full(len(chili), 3, dtype=np.int8)
    v = np.asarray(chili, dtype="float64")
    finite = np.isfinite(v)
    idx[finite & (v < CHILI_COOL_MAX)] = 0
    idx[finite & (v >= CHILI_COOL_MAX) & (v <= CHILI_WARM_MIN)] = 1
    idx[finite & (v > CHILI_WARM_MIN)] = 2
    return idx


def _bin_index(values, edges):
    """Index of the bin containing each value (-1 outside / nodata)."""
    v = np.asarray(values, dtype="float64")
    idx = np.searchsorted(edges, v, side="right") - 1
    idx[(v < edges[0]) | (v >= edges[-1]) | ~np.isfinite(v)] = -1
    return idx


def tile_partials(df, water_years, filter_tags=None, unit_types=None):
    """Partial sums for one tile's pixel table (``datacube.tabulate_tile``
    output: ints with -9999 nodata, floats with NaN), for every filter tag
    and unit type. Returns one long DataFrame; rows are keyed by
    (filter_tag, unit_type, unit_id, elevation, aspect, latitude,
    chili_class) with the geometric bins the unit type does not use left
    NaN. Reducible statistics only — see the module docstring.
    """
    water_years = [int(y) for y in water_years]
    filter_tags = list(FILTERS) if filter_tags is None else list(filter_tags)
    unit_types = list(UNIT_TYPES) if unit_types is None else list(unit_types)
    wy_cols = [f"runoff_onset_WY{y}" for y in water_years]

    # per-pixel derived columns, computed once for the whole tile (float32:
    # the values are small; np.bincount accumulates the sums in float64)
    base = pd.DataFrame(index=df.index)
    base["median"] = df["runoff_onset_median"].astype("float32")
    base["mad"] = df["runoff_onset_mad"].astype("float32")
    base["tres"] = df["temporal_resolution_median"].astype("float32")
    base["chili"] = df["chili"].astype("float32")
    fcf = df["forest_cover_fraction"].astype("float32")
    base["fcf"] = fcf.where(fcf >= 0)
    onset = df[wy_cols].astype("float32").where(df[wy_cols] > 0)
    base["nyears"] = onset.notna().sum(axis=1).astype("float32")
    for y, c in zip(water_years, wy_cols):
        base[f"on_{y}"] = onset[c]
        base[f"an_{y}"] = onset[c] - base["median"]
    del onset
    base["elevation"] = _bin_index(df["dem"].where(df["dem"] != NODATA), DEM_EDGES)
    aspect = df["aspect"].where(df["aspect"] != NODATA) % 360   # 360 -> 0 (north)
    base["aspect"] = _bin_index(aspect, ASPECT_EDGES)
    base["latitude"] = _bin_index(df["original_lat"], LAT_EDGES)
    base["chili_class"] = chili_class_index(base["chili"].values)

    out = []
    for tag in filter_tags:
        passing = np.zeros(len(df), dtype=bool)
        passing[df.index.get_indexer(apply_filter(df, tag).index)] = True
        for unit_type in unit_types:
            spec = UNIT_TYPES[unit_type]
            ids = df[spec["id_col"]].values
            ok = passing & (ids != NODATA)
            for b in spec["bins"]:
                if b not in spec["optional_bins"]:   # an optional bin may be -1 (undefined)
                    ok &= base[b].values >= 0
            if not ok.any():
                continue
            pos = np.flatnonzero(ok)   # row positions; columns are gathered one at a time
            keys = pd.DataFrame({"unit_id": ids[pos].astype("int64")})
            for b in spec["bins"]:
                keys[b] = base[b].values[pos]
            keys["chili_class"] = base["chili_class"].values[pos]
            part = _partial_sums(base, pos, keys, water_years)
            for b in ("elevation", "aspect", "latitude"):
                if b not in spec["bins"]:
                    part.insert(1, b, np.nan)
                else:  # bin index -> bin center; index -1 (undefined) -> NaN
                    centers = np.append(bin_centers(b), np.nan)
                    part[b] = centers[part[b].values.astype(int)]
            part.insert(0, "unit_type", unit_type)
            part.insert(0, "filter_tag", tag)
            out.append(part)
    if not out:
        return pd.DataFrame()
    result = pd.concat(out, ignore_index=True)
    result["chili_class"] = result["chili_class"].astype("int8")
    return result


def _partial_sums(base, pos, keys, water_years):
    """Reducible sums of the rows ``pos`` of ``base`` grouped by ``keys``.
    Memory-lean: keys are factorized once to integer group codes and every
    sum is an ``np.bincount`` over them; columns are gathered one at a time
    (the earlier pandas version built a 2M x 85 frame per filter x unit
    type and was OOM-killed on the dev box)."""
    kcols = list(keys.columns)
    codes, uniques = pd.factorize(pd.MultiIndex.from_frame(keys), sort=False)
    ngroups = len(uniques)
    out = pd.DataFrame(np.array([list(u) for u in uniques]), columns=kcols)
    for c in kcols:
        out[c] = out[c].astype("int64")

    def col(name):
        return base[name].values[pos].astype("float64")

    def bsum(values):
        v = np.asarray(values, dtype="float64")
        return np.bincount(codes, weights=np.where(np.isfinite(v), v, 0.0), minlength=ngroups)

    def bcount(values):
        return np.bincount(codes, weights=np.isfinite(np.asarray(values, dtype="float64")).astype("float64"),
                           minlength=ngroups)

    s = {"median": col("median"), "mad": col("mad"), "tres": col("tres"), "nyears": col("nyears")}
    med = s["median"]
    chili = col("chili")
    fcf = col("fcf")
    out["n"] = np.bincount(codes, minlength=ngroups).astype("float64")
    out["s_med"] = bsum(med); out["ss_med"] = bsum(med ** 2)
    out["n_mad"] = bcount(s["mad"]); out["s_mad"] = bsum(s["mad"]); out["ss_mad"] = bsum(s["mad"] ** 2)
    out["n_tres"] = bcount(s["tres"]); out["s_tres"] = bsum(s["tres"])
    out["s_nyears"] = bsum(s["nyears"])
    # correlation sums (pixel-level Pearson r of x vs the pixel median),
    # restricted to pixels where x is defined
    for x, name in ((chili, "chili"), (fcf, "fcf")):
        has = np.isfinite(x)
        out[f"n_{name}"] = bcount(x)
        out[f"s_{name}"] = bsum(x); out[f"ss_{name}"] = bsum(x ** 2)
        out[f"s_{name}_med"] = bsum(x * med)
        out[f"s_med_{name}"] = bsum(np.where(has, med, np.nan))
        out[f"ss_med_{name}"] = bsum(np.where(has, med ** 2, np.nan))
    for y in water_years:
        on = col(f"on_{y}"); an = col(f"an_{y}")
        out[f"n_{y}"] = bcount(on)
        out[f"s_{y}"] = bsum(on); out[f"ss_{y}"] = bsum(on ** 2)
        out[f"s_an_{y}"] = bsum(an); out[f"ss_an_{y}"] = bsum(an ** 2)
    return out


# ---------------------------------------------------------------------------
# REDUCE: partial sums (any number of tiles) -> one cube per unit type

def _mean_std(s, ss, n):
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s / n
        var = ss / n - mean ** 2
    std = np.sqrt(np.clip(var, 0, None))
    return mean, std


def _pearson(n, sx, sy, sxy, sxx, syy):
    with np.errstate(invalid="ignore", divide="ignore"):
        num = n * sxy - sx * sy
        den = np.sqrt((n * sxx - sx ** 2) * (n * syy - sy ** 2))
        r = num / den
    r[n < 3] = np.nan
    return r


def reduce_partials(partials, group, filter_tag, water_years):
    """Sum the partials of one output ``group`` (:data:`GROUPS`) and
    ``filter_tag`` over all tiles and build the cube: dims ``(<unit_dim>,
    <bins...>, chili_class, water_year)``; float32 means/stds, int32 counts;
    NaN where empty. Bins the group does not keep are summed out; basin ids
    are truncated to the group's HydroBASINS level (exact: the sums are
    additive and Pfafstetter codes nest by prefix); rows whose kept bin is
    undefined (flat pixels in ``continents_aspect``) are dropped. Continents
    are named (integer codes -> names, Australia -> Oceania).
    """
    gspec = GROUPS[group]
    unit_type = gspec["unit_type"]
    spec = dict(UNIT_TYPES[unit_type], bins=tuple(gspec["bins"]))
    if not set(spec["bins"]) <= set(UNIT_TYPES[unit_type]["bins"]):
        raise ValueError(f"group {group} keeps bins the map does not key: {spec['bins']}")
    water_years = [int(y) for y in water_years]
    sel = partials[(partials["unit_type"] == unit_type)
                   & (partials["filter_tag"] == filter_tag)]
    if sel.empty:
        raise ValueError(f"no partials for {unit_type} / {filter_tag}")
    sel = sel.dropna(subset=list(spec["bins"]))
    if "basin_level" in gspec:
        sel = sel.assign(unit_id=basin_level_ids(sel["unit_id"].values, gspec["basin_level"]))
    keys = ["unit_id", *spec["bins"], "chili_class"]
    num_cols = [c for c in sel.columns
                if c not in ("filter_tag", "unit_type", "elevation", "aspect",
                             "latitude", "chili_class", "unit_id", "tile_row", "tile_col")]
    sums = sel.groupby(keys, sort=True)[num_cols].sum(min_count=1).reset_index()

    off_axis = {}
    if unit_type == "continents":
        names = sums["unit_id"].map(CONTINENTS_ENUM).map(lambda n: CONTINENT_MERGE.get(n, n))
        # pixels outside the six study continents are dropped here, counted and recorded:
        # Antarctica (the dataset's southern edge nicks the peninsula) and any nodata code
        # (a tile without a continent polygon stored -9999, once wrapped to -15 in int8)
        off = ~names.isin(CONTINENT_NAMES)
        if off.any():
            off_axis = {str(CONTINENTS_ENUM.get(int(u), f"code {int(u)}")): int(n)
                        for u, n in sums.loc[off].groupby("unit_id")["n"].sum().items()}
            print(f"  continents: dropping pixels outside the six-continent axis: {off_axis}", flush=True)
            sums, names = sums.loc[~off], names.loc[~off]
        sums = sums.assign(unit_id=names.values)
        sums = sums.groupby(keys, sort=True)[num_cols].sum(min_count=1).reset_index()

    # the continent axis is the fixed six-continent set (empty ones stay NaN)
    # so panels indexed by continent name run on a partially processed version;
    # ranges and basins are whatever the partials contain
    units = (np.array(CONTINENT_NAMES) if unit_type == "continents"
             else np.array(sorted(sums["unit_id"].unique())))
    coords = {spec["unit_dim"]: units}
    for b in spec["bins"]:
        coords[b] = bin_centers(b)
    coords["chili_class"] = np.array(CHILI_CLASSES)
    coords["water_year"] = np.array(water_years)
    dims_static = (spec["unit_dim"], *spec["bins"], "chili_class")
    shape_static = tuple(len(coords[d]) for d in dims_static)

    idx = [pd.Index(units).get_indexer(sums["unit_id"])]
    for b in spec["bins"]:
        idx.append(pd.Index(coords[b]).get_indexer(sums[b]))
    idx.append(sums["chili_class"].values.astype(int))
    idx = tuple(idx)
    assert all((i >= 0).all() for i in idx), "partials bin or unit outside the coordinate grid"

    def dense(values, dtype=np.float32, fill=np.nan):
        arr = np.full(shape_static, fill, dtype=dtype)
        arr[idx] = values
        return arr

    n = sums["n"].values
    ds = xr.Dataset(coords=coords)
    mean, std = _mean_std(sums["s_med"].values, sums["ss_med"].values, n)
    ds["runoff_onset_median"] = (dims_static, dense(mean))
    ds["runoff_onset_median_std"] = (dims_static, dense(std))
    ds["runoff_onset_median_n"] = (dims_static, dense(n, np.int32, 0))
    mean, std = _mean_std(sums["s_mad"].values, sums["ss_mad"].values, sums["n_mad"].values)
    ds["runoff_onset_mad"] = (dims_static, dense(mean))
    ds["runoff_onset_mad_std"] = (dims_static, dense(std))
    ds["runoff_onset_mad_n"] = (dims_static, dense(sums["n_mad"].values, np.int32, 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        ds["temporal_resolution_median"] = (dims_static, dense(sums["s_tres"].values / sums["n_tres"].values))
        ds["n_years"] = (dims_static, dense(sums["s_nyears"].values / n))
    ds["temporal_resolution_median_n"] = (dims_static, dense(sums["n_tres"].values, np.int32, 0))

    dims_yearly = dims_static + ("water_year",)
    yr = {k: np.full(shape_static + (len(water_years),), np.nan, dtype=np.float32)
          for k in ("on", "on_std", "an", "an_std")}
    yr_n = np.zeros(shape_static + (len(water_years),), dtype=np.int32)
    for j, y in enumerate(water_years):
        ny = sums[f"n_{y}"].values
        m, sd = _mean_std(sums[f"s_{y}"].values, sums[f"ss_{y}"].values, ny)
        yr["on"][idx + (j,)] = m
        yr["on_std"][idx + (j,)] = sd
        m, sd = _mean_std(sums[f"s_an_{y}"].values, sums[f"ss_an_{y}"].values, ny)
        yr["an"][idx + (j,)] = m
        yr["an_std"][idx + (j,)] = sd
        yr_n[idx + (j,)] = ny
    ds["runoff_onset"] = (dims_yearly, yr["on"])
    ds["runoff_onset_std"] = (dims_yearly, yr["on_std"])
    ds["runoff_onset_n"] = (dims_yearly, yr_n)
    ds["runoff_onset_anomaly"] = (dims_yearly, yr["an"])
    ds["runoff_onset_anomaly_std"] = (dims_yearly, yr["an_std"])

    # pixel-level correlations of CHILI / FCF with the pixel median, per
    # geometric bin (summed over chili classes — the class is not a
    # dimension of a correlation)
    dims_corr = (spec["unit_dim"], *spec["bins"])
    csum = sums.groupby(["unit_id", *spec["bins"]], sort=True)[num_cols].sum(min_count=1).reset_index()
    cidx = [pd.Index(units).get_indexer(csum["unit_id"])]
    for b in spec["bins"]:
        cidx.append(pd.Index(coords[b]).get_indexer(csum[b]))
    cidx = tuple(cidx)
    shape_corr = tuple(len(coords[d]) for d in dims_corr)
    for x in ("chili", "fcf"):
        r = _pearson(csum[f"n_{x}"].values, csum[f"s_{x}"].values, csum[f"s_med_{x}"].values,
                     csum[f"s_{x}_med"].values, csum[f"ss_{x}"].values, csum[f"ss_med_{x}"].values)
        arr = np.full(shape_corr, np.nan, dtype=np.float32); arr[cidx] = r
        ds[f"{x}_corr"] = (dims_corr, arr)
        arr = np.zeros(shape_corr, dtype=np.int32); arr[cidx] = csum[f"n_{x}"].values
        ds[f"{x}_corr_n"] = (dims_corr, arr)

    for b in spec["bins"]:
        ds[b].attrs["units"] = BIN_UNITS[b]
        ds[b].attrs["bin_edges"] = BIN_EDGES[b].tolist()
    ds["water_year"].attrs["long_name"] = "water year"
    for v in ("runoff_onset_median", "runoff_onset"):
        ds[v].attrs["units"] = "day of water year"
    for v in ("runoff_onset_mad", "runoff_onset_anomaly", "runoff_onset_median_std",
              "runoff_onset_std", "runoff_onset_anomaly_std", "runoff_onset_mad_std"):
        ds[v].attrs["units"] = "days"
    ds["temporal_resolution_median"].attrs["units"] = "days"
    ds["n_years"].attrs["long_name"] = "mean number of valid water years per pixel"
    ds["chili_corr"].attrs["long_name"] = "Pearson r of CHILI vs pixel median onset, per bin"
    ds["fcf_corr"].attrs["long_name"] = "Pearson r of forest cover fraction vs pixel median onset, per bin"
    ds.attrs.update({
        "group": group, "unit_type": unit_type, "unit_id_column": spec["id_col"],
        "filter_tag": filter_tag, "filter_predicates": str(FILTERS[filter_tag]),
        "chili_class_thresholds": f"cool < {CHILI_COOL_MAX} <= neutral <= {CHILI_WARM_MIN} < warm; none = no CHILI",
        "statistics": "<var> = bin mean of pixel values; <var>_std = bin std; <var>_n = pixel count "
                      "(runoff_onset_anomaly shares runoff_onset_n). Marginals via aggregate.weighted_mean.",
    })
    if "basin_level" in gspec:
        ds.attrs["hydrobasins_level"] = gspec["basin_level"]
        ds.attrs["hydrobasins_level_stored"] = BASIN_LEVEL_STORED
        ds["river_basin"].attrs["long_name"] = f"HydroBASINS level-{gspec['basin_level']} PFAF_ID"
    if off_axis:
        ds.attrs["dropped_pixels_outside_continent_axis"] = str(off_axis)
    dropped = [b for b in UNIT_TYPES[unit_type]["bins"] if b not in spec["bins"]]
    if dropped:
        ds.attrs["summed_out_bins"] = ",".join(dropped)
    return ds


# ---------------------------------------------------------------------------
# read-side helpers (the notebooks' vocabulary for the unified schema)

# quantities whose counts live under another variable's name
_N_ALIAS = {"runoff_onset_anomaly": "runoff_onset_n",
            "runoff_onset_anomaly_std": "runoff_onset_n",
            "runoff_onset_std": "runoff_onset_n",
            "runoff_onset_median_std": "runoff_onset_median_n",
            "runoff_onset_mad_std": "runoff_onset_mad_n",
            "n_years": "runoff_onset_median_n"}


def counts_for(ds, var):
    """The pixel-count variable that weights ``var``."""
    return ds[_N_ALIAS.get(var, f"{var}_n")]


def weighted_mean(ds, var, dim):
    """Pixel-count-weighted mean of a bin-mean variable over ``dim`` (str or
    list) — the exact mean over all pixels in the merged bins. NaN where
    no pixels. Replaces the old ``basin_mean`` / ``runoff_onset_mean_anomaly``
    precomputes: e.g. ``weighted_mean(ds, 'runoff_onset_anomaly',
    ['elevation', 'aspect', 'chili_class'])``."""
    n = counts_for(ds, var).where(ds[var].notnull(), 0)
    return (ds[var] * n).sum(dim, min_count=1) / n.sum(dim).where(lambda x: x > 0)


def collapse(ds, dim="chili_class"):
    """Collapse a dimension exactly: means become count-weighted means,
    stds pooled, counts summed; variables without ``dim`` pass through.
    ``collapse(ds)`` gives the elevation x aspect (x latitude) cube the
    pre-2026 analyses worked with."""
    out = xr.Dataset(coords={k: v for k, v in ds.coords.items() if dim not in v.dims})
    for name, da in ds.data_vars.items():
        if dim not in da.dims:
            out[name] = da
            continue
        if name.endswith("_n") or name in ("runoff_onset_n",):
            out[name] = da.sum(dim)
        elif _N_ALIAS.get(name, f"{name}_n") not in ds:
            raise ValueError(f"cannot collapse {name!r}: no pixel count to weight it")
        elif name.endswith("_std"):
            base = name[:-4]
            n = counts_for(ds, base).where(ds[base].notnull(), 0)
            mean = weighted_mean(ds, base, dim)
            # pooled variance: sum n (sigma^2 + mu^2) / sum n - mu_bar^2
            pooled = ((n * (da ** 2 + ds[base] ** 2)).sum(dim, min_count=1)
                      / n.sum(dim).where(lambda x: x > 0) - mean ** 2)
            out[name] = np.sqrt(pooled.clip(min=0))
        else:
            out[name] = weighted_mean(ds, name, dim)
        out[name].attrs = da.attrs
    out.attrs = dict(ds.attrs, collapsed=dim)
    return out


def threshold(ds, min_n, variables=None):
    """Mask bin means (and stds) whose pixel count is below ``min_n`` —
    the old ``ds.where(ds.sel(statistic='count') > thresh)``."""
    out = ds.copy()
    for name, da in ds.data_vars.items():
        if variables is not None and name not in variables:
            continue
        if name.endswith("_n") or name in ("runoff_onset_n",) or name.endswith("_corr"):
            continue
        if _N_ALIAS.get(name, f"{name}_n") not in ds:   # e.g. the ERA5 zonal variables
            continue
        out[name] = da.where(counts_for(ds, name) >= min_n)
    return out


def elevation_relative(da, dim="aspect"):
    """Deviation from the per-elevation median across ``dim`` (the triplets'
    third panel; was the stored ``runoff_onset_elev_relative``)."""
    return da - da.median(dim)


# ---------------------------------------------------------------------------
# the partials cache: the fleet product, downloaded once per tile

def sync_partials(config, cache_dir, workers=16, log=print):
    """Download the partials parquets of ``config.version`` that are missing from ``cache_dir``
    (``paths.partials_cache(version)``) and delete cached files Azure no longer lists (a redone
    tile). Re-runs only fetch new tiles. Needs the Azure SAS token. Returns the sorted local paths."""
    import os
    from concurrent.futures import ThreadPoolExecutor
    fs = config.azure_blob_fs
    prefix = f"{settings.PARTIALS_PREFIX}/{config.version}"
    remote = {p.rsplit('/', 1)[-1]: p for p in fs.ls(prefix, detail=False) if p.endswith('.parquet')}
    local = {f for f in os.listdir(cache_dir) if f.endswith('.parquet')}
    stale = local - set(remote)
    for f in stale:
        os.remove(cache_dir / f)
    missing = sorted(set(remote) - local)
    log(f"partials on Azure: {len(remote)} | cached: {len(local) - len(stale)} | "
        f"downloading {len(missing)} | dropped stale {len(stale)}")

    def fetch(name):
        fs.get_file(remote[name], str(cache_dir / f"{name}.part"))
        os.replace(cache_dir / f"{name}.part", cache_dir / name)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(fetch, missing))
    return sorted(cache_dir / f for f in remote)
