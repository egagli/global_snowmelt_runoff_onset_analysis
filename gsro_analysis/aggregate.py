"""The aggregation engine: per-tile pixel tables -> per-tile PARTIAL SUMS
(the fleet "map", :func:`tile_partials`) -> one SPARSE TABLE per unit type
(the "reduce", :func:`reduce_sparse`, DuckDB over the partials parquets) ->
dense cubes on demand (:func:`cube` / :func:`view`), plus the read-side
helpers every analysis notebook uses (:func:`weighted_mean`,
:func:`collapse`, :func:`threshold`, :func:`elevation_relative`).

Design (partials schema 3, 2026-10-20; schema 1/2 history in docs/aggregation_lineage.md):

- Every statistic is *reducible*: means and standard deviations from
  (n, sum, sum of squares), correlations from the six Pearson sums, onset
  distributions from fixed-bin histogram counts. Partials sum exactly over
  tiles, so the reduce is a GROUP BY ... SUM.
- **Which pixels.** The pixel table (``datacube.tabulate_tile``) holds every
  pixel of the tile's 80 m UTM grid with a DEM (~6,400 m2 each: counts are
  area). Two predicate sets act on it in the map:
  ``TERRAIN`` (land: not ESA WorldCover water) gives ``n_terrain``, the
  terrain hypsometry; ``ANALYSED`` (a valid 11-yr median onset, not the
  Liston & Sturm "ephemeral" snow class, not built-up or water) gives every
  onset statistic. Forest cover is no longer a filter: it is a KEY.
- **Keys of a row:** ``unit_type``, ``unit_id``, the unit type's geometric
  bins (elevation 100 m, aspect 15 deg — NaN for flat pixels —, latitude 1 deg
  for continents) and three class keys on every unit type:
  ``fcf_class`` (forest cover 0-10 %, (10, 20], ..., (90, 100], none = no FCF
  value: :data:`FCF_CLASSES`), ``slope_class`` (0-10, ..., 40-50, 50+ deg:
  :data:`SLOPE_CLASSES`) and ``perennial`` (1 = ESA WorldCover class 70,
  permanent snow and ice). CHILI is no longer a key (its sums stay, so the mean
  CHILI and its correlation with onset are still available per bin).
  The two former filter tags are views: ``fcf_lte_50`` = FCF classes 0-4,
  ``full_dataset`` = every class including "none" (:data:`FILTER_VIEWS`).
- **Histogram rows** (``onset_bin`` not null): for :data:`HIST_UNIT_TYPES`,
  per unit x elevation x perennial x 5-day onset bin, counts of the 11-yr
  median and of each year's onset over analysed pixels with FCF <= 50 %.
- **Water towers** (Immerzeel et al. 2020): unit ids are signed by the
  ancillary layer ``WTU_SIGNED`` (the WTU id inside the tower polygon, id +
  1000 in its downstream part), so towers and whole basins come from one map.
- **Reduce.** :func:`reduce_sparse` sums one unit type over all tiles with
  DuckDB (out of core: ~17 GB of partials per version) into a sparse table,
  one row per non-empty key combination, written straight to parquet (the
  continents table alone has ~7 M rows x 79 sums, too large for pandas on a
  16 GB machine). :func:`load_sparse` reads it back summed down to the keys a
  view needs. :func:`view` turns a sparse table
  into the dense cube a notebook needs, keeping the keys it names and summing
  the others; :func:`cube` builds the named default cubes (:data:`GROUPS`)
  with the old variable names (no ``chili_class`` dimension any more).
- ``water_years`` is always a parameter (``config.water_years``).
- No cluster service, no dask: the map runs inside the fleet worker, the
  reduce is DuckDB + pandas + numpy in a notebook.
"""

import math
import operator
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from gsro_analysis import settings

# ---------------------------------------------------------------------------
# pixel predicates (AND-ed on the pixel table)
TERRAIN = [("esa_worldcover", "!=", 80)]                      # land: every pixel with a DEM that is not water
ANALYSED = [("runoff_onset_median", ">", 0),                  # a valid 11-yr median onset (-9999 = none)
            ("snow_classification", "!=", 4),                 # not Liston & Sturm 'ephemeral (includes no snow)'
            ("esa_worldcover", "not in", [50, 80])]           # not built-up, not water

_OPS = {
    "==": operator.eq, "!=": operator.ne,
    ">": operator.gt, ">=": operator.ge, "<": operator.lt, "<=": operator.le,
    "in": lambda s, v: s.isin(v), "not in": lambda s, v: ~s.isin(v),
}


def predicate_mask(df, predicates):
    """Boolean mask of the rows of a pixel table passing AND-ed predicates."""
    mask = np.ones(len(df), dtype=bool)
    for col, op, val in predicates:
        mask &= np.asarray(_OPS[op](df[col], val))
    return mask


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

# CHILI insolation classes (Theobald et al. 2015 thresholds): no longer a partials key (schema 3);
# kept for chili_class_index on pixel tables
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

# class keys (schema 3): integer codes in the partials, labels as cube coordinates
FCF_CLASSES = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "none"]
SLOPE_CLASSES = ["0-10", "10-20", "20-30", "30-40", "40-50", "50+"]
PERENNIAL_CLASSES = ["seasonal", "perennial"]
CLASS_KEYS = {"fcf_class": FCF_CLASSES, "slope_class": SLOPE_CLASSES, "perennial": PERENNIAL_CLASSES}
FCF_NONE = FCF_CLASSES.index("none")
PERENNIAL_WORLDCOVER_CLASS = 70        # ESA WorldCover 2021 'permanent snow and ice'
# the former filter tags, now selections on fcf_class
FILTER_VIEWS = {"fcf_lte_50": list(range(5)), "full_dataset": list(range(len(FCF_CLASSES)))}
FILTERS = FILTER_VIEWS                  # the name the notebooks iterate over
# onset histograms: 5-day bins of the day of water year
ONSET_BIN_EDGES = np.arange(0, 366 + 5, 5)
HIST_UNIT_TYPES = ("mountain_ranges", "grdc_basins", "water_towers")
HIST_MAX_FCF_CLASS = 4                  # histograms count analysed pixels with FCF <= 50 %


def bin_centers(name):
    e = BIN_EDGES[name]
    return e[:-1] + np.diff(e) / 2


# The unit types the MAP keys partial rows by (geometric bins; aspect is NaN for flat pixels and only the
# views that keep aspect drop them). ``id_col`` is the ancillary layer: GMBA_V2_ID (standard 300 inventory),
# PFAF_ID at HydroBASINS level :data:`BASIN_LEVEL_STORED`, the USGS continent code, GRDC_MRBID (GRDC Major
# River Basins 2020), WTU_SIGNED (Immerzeel et al. 2020 water towers: id inside, id + 1000 downstream).
UNIT_TYPES = {
    "mountain_ranges": {"id_col": "GMBA_V2_ID", "unit_dim": "mountain_range", "bins": ("elevation", "aspect")},
    "river_basins":    {"id_col": "PFAF_ID", "unit_dim": "river_basin", "bins": ("elevation",)},
    "continents":      {"id_col": "continent", "unit_dim": "continent", "bins": ("latitude", "elevation", "aspect")},
    "grdc_basins":     {"id_col": "GRDC_MRBID", "unit_dim": "grdc_basin", "bins": ("elevation", "aspect")},
    "water_towers":    {"id_col": "WTU_SIGNED", "unit_dim": "water_tower", "bins": ("elevation", "aspect")},
}
OPTIONAL_BINS = ("aspect",)
BASIN_LEVEL_STORED = settings.BASIN_ATLAS_STORED_LEVEL   # 6
WTU_DOWNSTREAM_OFFSET = 1000

# The default cubes (``cube(sparse, group, filter_tag, ...)``): a unit type, the bins kept (the others and
# the class keys are summed out — exact), and for basins the HydroBASINS level (a digit prefix of the stored
# id). The aggregation notebooks name the file ``all_<group>_<filter_tag>.nc``. Any other combination
# (e.g. slope_class kept) is one ``view()`` call on the sparse table.
GROUPS = {
    "mountain_ranges":   {"unit_type": "mountain_ranges", "bins": ("elevation", "aspect")},
    "river_basins":      {"unit_type": "river_basins", "bins": ("elevation",), "basin_level": 5},
    "river_basins_l6":   {"unit_type": "river_basins", "bins": ("elevation",), "basin_level": 6},
    "continents":        {"unit_type": "continents", "bins": ("latitude", "elevation")},
    "continents_aspect": {"unit_type": "continents", "bins": ("latitude", "elevation", "aspect")},
    "grdc_basins":       {"unit_type": "grdc_basins", "bins": ("elevation",)},
    "grdc_basins_aspect": {"unit_type": "grdc_basins", "bins": ("elevation", "aspect")},
    "water_towers":      {"unit_type": "water_towers", "bins": ("elevation", "aspect")},
}
DEFAULT_GROUPS = ("mountain_ranges", "river_basins", "continents", "grdc_basins", "water_towers")

# Version of the partials layout, stamped into every partials parquet's key-value metadata
# (``gsro_partials_schema``). The fleet re-maps a tile whose blob carries an older stamp (or none):
# 1 = 2026-09-03 (three unit types), 2 = 2026-10-09 (+ grdc_basins), 3 = 2026-10-20 (terrain counts, FCF /
# slope / perennial keys replace the filter tags and CHILI, latitude moments, water towers, onset histograms).
PARTIALS_SCHEMA = 3

KEY_COLS = ["unit_type", "unit_id", "elevation", "aspect", "latitude", "fcf_class", "slope_class", "perennial",
            "onset_bin"]


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


def fcf_class_index(fcf):
    """Forest-cover class per :data:`FCF_CLASSES`: 0 = [0, 10], k = (10k, 10k + 10] for k = 1..9 (so the
    former ``fcf_lte_50`` filter, 0 <= FCF <= 50, is exactly classes 0-4), 10 = no FCF value (< 0 / NaN)."""
    v = np.asarray(fcf, dtype="float64")
    idx = np.full(len(v), FCF_NONE, dtype=np.int8)
    ok = np.isfinite(v) & (v >= 0)
    idx[ok] = np.clip(np.where(v[ok] <= 10, 0, np.ceil(v[ok] / 10) - 1), 0, 9)
    return idx


def slope_class_index(slope):
    """Slope class per :data:`SLOPE_CLASSES` (10 deg classes, 50+ the last; a missing slope counts as 0-10)."""
    v = np.asarray(slope, dtype="float64")
    return np.where(np.isfinite(v) & (v >= 0), np.clip(np.floor(v / 10), 0, 5), 0).astype(np.int8)


def _bin_index(values, edges):
    """Index of the bin containing each value (-1 outside / nodata)."""
    v = np.asarray(values, dtype="float64")
    idx = np.searchsorted(edges, v, side="right") - 1
    idx[(v < edges[0]) | (v >= edges[-1]) | ~np.isfinite(v)] = -1
    return idx


def _num(df, col):
    """A pixel-table column as float32 with the -9999 fill as NaN."""
    v = df[col].to_numpy(dtype="float32", copy=True)
    v[v == NODATA] = np.nan
    return v


def tile_partials(df, water_years, unit_types=None):
    """Partial sums for one tile's pixel table (``datacube.tabulate_tile``:
    every pixel with a DEM; ints with -9999 nodata, floats with NaN) for every
    unit type. Returns one long DataFrame keyed by :data:`KEY_COLS` (bins the
    unit type does not use are NaN; ``onset_bin`` is NaN on sum rows and set
    on histogram rows). Reducible statistics only — see the module docstring.
    """
    water_years = [int(y) for y in water_years]
    unit_types = list(UNIT_TYPES) if unit_types is None else list(unit_types)
    wy_cols = [f"runoff_onset_WY{y}" for y in water_years]
    land = predicate_mask(df, TERRAIN)
    analysed = predicate_mask(df, ANALYSED)

    # per-pixel derived columns, computed once for the whole tile; onset statistics only on analysed pixels
    base = pd.DataFrame(index=df.index)
    med = _num(df, "runoff_onset_median"); med[~analysed] = np.nan
    base["median"] = med
    base["mad"] = np.where(analysed, df["runoff_onset_mad"].to_numpy("float32"), np.nan)
    base["tres"] = np.where(analysed, df["temporal_resolution_median"].to_numpy("float32"), np.nan)
    base["chili"] = np.where(analysed, df["chili"].to_numpy("float32"), np.nan)
    fcf = _num(df, "forest_cover_fraction"); fcf[fcf < 0] = np.nan
    base["fcf"] = np.where(analysed, fcf, np.nan)
    base["lat"] = np.where(analysed, df["original_lat"].to_numpy("float64"), np.nan)
    onset = df[wy_cols].to_numpy("float32", copy=True)
    onset[(onset <= 0) | ~analysed[:, None]] = np.nan
    base["nyears"] = np.where(analysed, np.isfinite(onset).sum(axis=1), np.nan).astype("float32")
    for j, y in enumerate(water_years):
        base[f"on_{y}"] = onset[:, j]
        base[f"an_{y}"] = onset[:, j] - med
    base["elevation"] = _bin_index(_num(df, "dem"), DEM_EDGES)
    aspect = _num(df, "aspect") % 360                           # 360 -> 0 (north); flat pixels NaN
    base["aspect"] = _bin_index(aspect, ASPECT_EDGES)
    base["latitude"] = _bin_index(df["original_lat"], LAT_EDGES)
    base["fcf_class"] = fcf_class_index(fcf)
    base["slope_class"] = slope_class_index(_num(df, "slope"))
    base["perennial"] = (df["esa_worldcover"].to_numpy() == PERENNIAL_WORLDCOVER_CLASS).astype(np.int8)
    onset_bin_med = _bin_index(med, ONSET_BIN_EDGES)
    onset_bin_yr = {y: _bin_index(onset[:, j], ONSET_BIN_EDGES) for j, y in enumerate(water_years)}

    out = []
    for unit_type in unit_types:
        spec = UNIT_TYPES[unit_type]
        ids = df[spec["id_col"]].to_numpy()
        ok = (ids != NODATA) & land
        for b in spec["bins"]:
            if b not in OPTIONAL_BINS:
                ok &= base[b].to_numpy() >= 0
        if not ok.any():
            continue
        pos = np.flatnonzero(ok)
        key_names = list(spec["bins"]) + list(CLASS_KEYS)
        keys = pd.DataFrame({"unit_id": ids[pos].astype("int64")})
        for k in key_names:
            keys[k] = base[k].to_numpy()[pos]
        part = _partial_sums(base, pos, keys, water_years, analysed[pos])
        for b in ("elevation", "aspect", "latitude"):
            if b not in spec["bins"]:
                part[b] = np.nan
            else:                       # bin index -> bin center; index -1 (undefined) -> NaN
                centers = np.append(bin_centers(b), np.nan)
                part[b] = centers[part[b].to_numpy().astype(int)]
        part["onset_bin"] = np.nan
        part.insert(0, "unit_type", unit_type)
        out.append(part)
        if unit_type in HIST_UNIT_TYPES:
            h = analysed[pos] & (base["fcf_class"].to_numpy()[pos] <= HIST_MAX_FCF_CLASS)
            if h.any():
                out.append(_histogram_rows(unit_type, ids[pos][h], base["elevation"].to_numpy()[pos][h],
                                           base["perennial"].to_numpy()[pos][h], onset_bin_med[pos][h],
                                           {y: v[pos][h] for y, v in onset_bin_yr.items()}))
    if not out:
        return pd.DataFrame()
    result = pd.concat(out, ignore_index=True)
    for k in CLASS_KEYS:
        result[k] = result[k].astype("int8")
    return result[KEY_COLS + [c for c in result.columns if c not in KEY_COLS]]


def _histogram_rows(unit_type, ids, elev_idx, perennial, med_bin, year_bins):
    """Onset histogram rows: counts per (unit, elevation, perennial, 5-day bin) of the 11-yr median
    (``h_med``) and of each year's onset (``h_<year>``)."""
    frames = []
    for name, b in [("h_med", med_bin)] + [(f"h_{y}", v) for y, v in year_bins.items()]:
        ok = b >= 0
        if not ok.any():
            continue
        frames.append(pd.DataFrame({"unit_id": ids[ok].astype("int64"), "elevation": elev_idx[ok],
                                    "perennial": perennial[ok], "onset_bin": b[ok]})
                      .value_counts().rename(name))
    h = pd.concat(frames, axis=1).fillna(0).reset_index()
    h["elevation"] = bin_centers("elevation")[h["elevation"].to_numpy().astype(int)]
    h["onset_bin"] = ONSET_BIN_EDGES[h["onset_bin"].to_numpy().astype(int)].astype("float64")   # bin start day
    h.insert(0, "unit_type", unit_type)
    for k in ("aspect", "latitude"):
        h[k] = np.nan
    h["fcf_class"] = -1                 # histograms are over FCF classes 0-4 together
    h["slope_class"] = -1               # and over all slopes
    return h


def _partial_sums(base, pos, keys, water_years, analysed):
    """Reducible sums of the rows ``pos`` of ``base`` grouped by ``keys``: ``n_terrain`` over all of them,
    the onset statistics over the ``analysed`` subset. Memory-lean: keys are factorized once to integer
    group codes and every sum is an ``np.bincount`` over them."""
    kcols = list(keys.columns)
    codes, uniques = pd.factorize(pd.MultiIndex.from_frame(keys), sort=False)
    ngroups = len(uniques)
    out = pd.DataFrame(np.array([list(u) for u in uniques]), columns=kcols)
    for c in kcols:
        out[c] = out[c].astype("int64")
    a = np.asarray(analysed)
    ca = codes[a]

    def col(name):
        return base[name].to_numpy()[pos][a].astype("float64")

    def bsum(values):
        v = np.asarray(values, dtype="float64")
        return np.bincount(ca, weights=np.where(np.isfinite(v), v, 0.0), minlength=ngroups)

    def bcount(values):
        return np.bincount(ca, weights=np.isfinite(np.asarray(values, dtype="float64")).astype("float64"),
                           minlength=ngroups)

    out["n_terrain"] = np.bincount(codes, minlength=ngroups).astype("float64")
    med = col("median")
    out["n"] = np.bincount(ca, minlength=ngroups).astype("float64")
    out["s_med"] = bsum(med); out["ss_med"] = bsum(med ** 2)
    mad, tres = col("mad"), col("tres")
    out["n_mad"] = bcount(mad); out["s_mad"] = bsum(mad); out["ss_mad"] = bsum(mad ** 2)
    out["n_tres"] = bcount(tres); out["s_tres"] = bsum(tres)
    out["s_nyears"] = bsum(col("nyears"))
    lat = col("lat")
    out["s_lat"] = bsum(lat); out["ss_lat"] = bsum(lat ** 2)
    # correlation sums (pixel-level Pearson r of x vs the pixel median), restricted to pixels where x is defined
    for name in ("chili", "fcf"):
        x = col(name)
        has = np.isfinite(x)
        out[f"n_{name}"] = bcount(x)
        out[f"s_{name}"] = bsum(x); out[f"ss_{name}"] = bsum(x ** 2)
        out[f"s_{name}_med"] = bsum(x * med)
        out[f"s_med_{name}"] = bsum(np.where(has, med, np.nan))
        out[f"ss_med_{name}"] = bsum(np.where(has, med ** 2, np.nan))
    for y in water_years:
        on, an = col(f"on_{y}"), col(f"an_{y}")
        out[f"n_{y}"] = bcount(on)
        out[f"s_{y}"] = bsum(on); out[f"ss_{y}"] = bsum(on ** 2)
        out[f"s_an_{y}"] = bsum(an); out[f"ss_an_{y}"] = bsum(an ** 2)
    return out


# ---------------------------------------------------------------------------
# REDUCE: partials of all tiles -> one sparse table per unit type (DuckDB) -> dense cubes

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


def _sum_columns(columns):
    return [c for c in columns if c not in KEY_COLS and c not in ("tile_row", "tile_col")]


def _schema3_files(partial_files):
    """The partials files to reduce: verified-empty tiles (no columns) are skipped; a file stamped with an
    older ``gsro_partials_schema`` raises (re-map the tile; mixing schemas would sum different keys)."""
    import pyarrow.parquet as pq
    keep = []
    for f in partial_files:
        sch = pq.read_schema(f)
        if "unit_type" not in sch.names:
            continue
        stamp = (sch.metadata or {}).get(b"gsro_partials_schema", b"?").decode()
        if stamp != str(PARTIALS_SCHEMA):
            raise ValueError(f"{f}: partials schema {stamp}, this code reduces schema {PARTIALS_SCHEMA}: re-map the tile")
        keep.append(str(f))
    return keep


# DuckDB's memory cap (its default is 80 % of RAM, which the reduce exceeds on a 16 GB machine) and the
# target rows per GROUP BY partition of :func:`reduce_sparse` (each partition's hash table stays ~1 GB)
DUCKDB_MEMORY_LIMIT = os.environ.get("GSRO_DUCKDB_MEMORY_LIMIT", "2GB")
DUCKDB_THREADS = int(os.environ.get("GSRO_DUCKDB_THREADS", "4"))   # buffers outside the memory limit grow per thread
SPARSE_ROWS_PER_PARTITION = 1_000_000
# files per first-stage batch: DuckDB holds the footers of every file it scans (~0.5 MB each), so the reduce
# sums batches of tiles into intermediate parquets first, then sums those
FILES_PER_BATCH = 500


def _duckdb_connect(threads=None):
    """An in-memory DuckDB connection with :data:`DUCKDB_MEMORY_LIMIT` and its own spill directory
    (returned; the caller removes it)."""
    import duckdb
    tmp = tempfile.mkdtemp(prefix="gsro_duckdb_")
    con = duckdb.connect()
    con.execute(f"SET memory_limit='{DUCKDB_MEMORY_LIMIT}'")
    con.execute(f"SET temp_directory='{tmp}'")
    con.execute("SET preserve_insertion_order=false")
    con.execute(f"SET threads={int(threads or DUCKDB_THREADS)}")
    return con, tmp


def reduce_sparse(partial_files, unit_type, histograms=False, threads=None, out_path=None):
    """Sum the partials of ``unit_type`` over every tile (DuckDB, out of core) into its SPARSE TABLE: one
    row per non-empty key combination (:data:`KEY_COLS`), the partial-sum columns summed. With
    ``histograms`` the onset-histogram rows instead (``onset_bin`` set). Empty (verified-empty) tiles and
    older-schema files without the unit type contribute nothing. The GROUP BY runs in partitions of unit
    ids of about :data:`SPARSE_ROWS_PER_PARTITION` rows each. With ``out_path`` the table is written to
    that parquet without passing through pandas and the path is returned (read it back with
    :func:`load_sparse`); without, a pandas DataFrame is returned (small inputs: tests, fixtures)."""
    files = _schema3_files(partial_files)
    keys = [k for k in KEY_COLS if k != "unit_type"]
    where = f"unit_type = '{unit_type}' AND onset_bin IS {'NOT ' if histograms else ''}NULL"
    con, tmp = _duckdb_connect(threads)

    def grouped(sources, where, sum_cols, prefix):
        """GROUP BY the keys over parquet ``sources`` in unit-id partitions; returns the partition files."""
        rel = con.read_parquet(sources, union_by_name=True)
        sum_cols = [c for c in sum_cols if c in rel.columns]   # a batch may lack some (e.g. histogram) columns
        n_keys = con.execute(f"SELECT approx_count_distinct(hash({', '.join(keys)})) FROM rel WHERE {where}").fetchone()[0]
        n_parts = max(1, math.ceil(n_keys / SPARSE_ROWS_PER_PARTITION))
        out = []
        for p in range(n_parts):
            part = f"((CAST(unit_id AS BIGINT) % {n_parts}) + {n_parts}) % {n_parts} = {p}"
            q = (f"SELECT '{unit_type}' AS unit_type, {', '.join(keys)}, "
                 f"{', '.join(f'SUM({c}) AS {c}' for c in sum_cols)} "
                 f"FROM rel WHERE {where} AND {part} GROUP BY {', '.join(keys)}")
            out.append(os.path.join(tmp, f"{prefix}_{p}.parquet"))
            con.execute(f"COPY ({q}) TO '{out[-1]}' (FORMAT parquet)")
        return out

    try:
        import pyarrow.parquet as pq
        cols = list(dict.fromkeys(c for f in files for c in pq.read_schema(f).names))   # the union of the schemas
        if "unit_type" not in cols:
            empty = pd.DataFrame(columns=KEY_COLS)
            if out_path is None:
                return empty
            empty.to_parquet(out_path, index=False)
            return Path(out_path)
        hist_cols = [c for c in cols if c.startswith("h_")]
        sum_cols = [c for c in _sum_columns(cols) if (c in hist_cols) == histograms]
        # stage 1: batches of tiles -> intermediate sums (only when there are more files than one batch)
        batches = [files[i:i + FILES_PER_BATCH] for i in range(0, len(files), FILES_PER_BATCH)]
        if len(batches) > 1:
            sources = [f for b, batch in enumerate(batches) for f in grouped(batch, where, sum_cols, f"batch{b}")]
            where = "TRUE"                       # the intermediates hold only the selected rows
        else:
            sources = files
        # stage 2: the final sums
        part_files = grouped(sources, where, sum_cols, "part")
        union = f"SELECT * FROM read_parquet({part_files!r}, union_by_name = true)"
        if out_path is None:
            return con.execute(union).df()
        con.execute(f"COPY ({union}) TO '{out_path}' (FORMAT parquet)")
        return Path(out_path)
    finally:
        con.close()
        shutil.rmtree(tmp, ignore_errors=True)


def load_sparse(path, keys, fcf=None, slope=None, perennial=None, threads=None):
    """A sparse table written by :func:`reduce_sparse` (``out_path``), summed down to ``unit_id`` plus
    ``keys`` (geometric bins and/or class keys) over the selection ``fcf`` / ``slope`` / ``perennial``
    (lists of codes; ``None`` = all), in DuckDB. The class keys summed out are recorded in
    ``df.attrs['presummed_selection']``, which :func:`view` checks against the selection it is asked
    for. Use it instead of reading the whole table: :func:`view` gives the same cube from either."""
    keys = list(dict.fromkeys(["unit_id", *keys]))
    con, tmp = _duckdb_connect(threads)
    try:
        rel = con.read_parquet(str(path))
        sum_cols = [c for c in _sum_columns(rel.columns) if not c.startswith("h_")]
        where, presummed = [], {}
        for key, codes in (("fcf_class", fcf), ("slope_class", slope), ("perennial", perennial)):
            if codes is not None:
                where.append(f"{key} IN ({', '.join(str(int(c)) for c in codes)})")
            if key not in keys:
                presummed[key] = None if codes is None else sorted(int(c) for c in codes)
        if sum_cols:
            q = (f"SELECT {', '.join(keys)}, {', '.join(f'SUM({c}) AS {c}' for c in sum_cols)} FROM rel "
                 f"{'WHERE ' + ' AND '.join(where) if where else ''} GROUP BY {', '.join(keys)}")
            out = con.execute(q).df()
        else:                                   # an empty table (no rows of the unit type)
            out = pd.DataFrame(columns=keys)
    finally:
        con.close()
        shutil.rmtree(tmp, ignore_errors=True)
    out.attrs["presummed_selection"] = presummed
    return out


def cube_from_sparse(path, group, filter_tag, water_years, **kw):
    """:func:`cube` from a sparse table on disk: :func:`load_sparse` with the group's bins kept and the
    filter tag's FCF selection applied."""
    g = GROUPS[group]
    sparse = load_sparse(path, g["bins"], fcf=FILTER_VIEWS[filter_tag])
    return cube(sparse, group, filter_tag, water_years, **kw)


def sparse_from_frame(partials, unit_type, histograms=False):
    """:func:`reduce_sparse` for partials already in memory (a DataFrame of one or more tiles; tests, the
    CI fixtures)."""
    sel = partials[(partials["unit_type"] == unit_type) & (partials["onset_bin"].notna() == histograms)]
    cols = [c for c in _sum_columns(sel.columns) if c.startswith("h_") == histograms]
    keys = [k for k in KEY_COLS if k != "unit_type"]
    out = sel.groupby(keys, dropna=False, sort=True)[cols].sum(min_count=1).reset_index()
    out.insert(0, "unit_type", unit_type)
    return out


def _unit_ids(sparse, unit_type, basin_level=None):
    """Unit ids of a sparse table as the cube's unit axis values: HydroBASINS ids truncated to
    ``basin_level``; continents named (Australia -> Oceania; codes outside the six-continent axis
    returned as None)."""
    ids = sparse["unit_id"].to_numpy().astype("int64")
    if unit_type == "river_basins" and basin_level:
        return basin_level_ids(ids, basin_level)
    if unit_type == "continents":
        names = pd.Series(ids).map(CONTINENTS_ENUM).map(lambda n: CONTINENT_MERGE.get(n, n))
        return names.where(names.isin(CONTINENT_NAMES), None).to_numpy()
    return ids


def view(sparse, unit_type, bins, water_years, fcf=FILTER_VIEWS["fcf_lte_50"], slope=None, perennial=None,
         keep=(), basin_level=None, analysed_units_only=True):
    """A dense cube from a sparse table (:func:`reduce_sparse`): dims ``(<unit_dim>, <bins...>, <keep...>,
    water_year)`` where ``keep`` names class keys to keep as dimensions (``'fcf_class'``,
    ``'slope_class'``, ``'perennial'``); every other key is summed out over the selection
    ``fcf`` (list of FCF class codes; default the analyses' FCF <= 50 %), ``slope`` and ``perennial``
    (lists of codes; ``None`` = all). Rows whose kept geometric bin is undefined (flat pixels when aspect is
    kept) are dropped, as in schema 2. Variables as the schema-2 cubes (``runoff_onset_median``,
    ``..._std``, ``..._n``, ``runoff_onset_mad``, ``temporal_resolution_median``, ``n_years``,
    ``runoff_onset``, ``runoff_onset_anomaly``, ``runoff_onset_n``, ``chili_corr``, ``fcf_corr``) plus
    ``terrain_n`` (land pixels of the same selection, onset or not), ``chili_mean`` and ``latitude_mean``.
    With ``analysed_units_only`` (default) the unit axis holds the units with analysed pixels in the
    selection, as in schema 2; ``False`` keeps snow-free units too (terrain only, every onset variable NaN)."""
    spec = UNIT_TYPES[unit_type]
    water_years = [int(y) for y in water_years]
    presummed = sparse.attrs.get("presummed_selection", {})     # class keys load_sparse already summed out
    sel = sparse
    for key, codes in (("fcf_class", fcf), ("slope_class", slope), ("perennial", perennial)):
        if key in presummed:
            if (None if codes is None else sorted(int(c) for c in codes)) != presummed[key]:
                raise ValueError(f"{key} was summed out over {presummed[key]} by load_sparse; this view selects {codes}")
            continue
        if codes is not None:
            sel = sel[sel[key].isin(list(codes))]
    sel = sel.dropna(subset=list(bins))
    units = _unit_ids(sel, unit_type, basin_level)
    sel = sel.assign(unit_id=units)
    sel = sel[pd.notna(sel["unit_id"])]
    if sel.empty:
        raise ValueError(f"no {unit_type} rows in the selection")
    keys = ["unit_id", *bins, *keep]
    num_cols = [c for c in _sum_columns(sel.columns) if not c.startswith("h_")]
    sums = sel.groupby(keys, sort=True)[num_cols].sum(min_count=1).reset_index()
    if analysed_units_only and unit_type != "continents":
        analysed = sums.groupby("unit_id")["n"].sum()
        sums = sums[sums["unit_id"].isin(analysed.index[analysed > 0])]
        if sums.empty:
            raise ValueError(f"no {unit_type} units with analysed pixels in the selection")

    unit_axis = (np.array(CONTINENT_NAMES) if unit_type == "continents" else np.array(sorted(sums["unit_id"].unique())))
    coords = {spec["unit_dim"]: unit_axis}
    for b in bins:
        coords[b] = bin_centers(b)
    for k in keep:
        coords[k] = np.array(CLASS_KEYS[k])
    coords["water_year"] = np.array(water_years)
    dims_static = (spec["unit_dim"], *bins, *keep)
    shape_static = tuple(len(coords[d]) for d in dims_static)
    idx = [pd.Index(unit_axis).get_indexer(sums["unit_id"])]
    for b in bins:
        idx.append(pd.Index(coords[b]).get_indexer(sums[b]))
    for k in keep:
        idx.append(sums[k].to_numpy().astype(int))
    idx = tuple(idx)
    assert all((i >= 0).all() for i in idx), "sparse-table bin or unit outside the coordinate grid"

    def dense(values, dtype=np.float32, fill=np.nan):
        arr = np.full(shape_static, fill, dtype=dtype)
        arr[idx] = values
        return arr

    n = sums["n"].to_numpy()
    nz = np.where(n > 0, n, np.nan)
    ds = xr.Dataset(coords=coords)
    mean, std = _mean_std(sums["s_med"].to_numpy(), sums["ss_med"].to_numpy(), nz)
    ds["runoff_onset_median"] = (dims_static, dense(mean))
    ds["runoff_onset_median_std"] = (dims_static, dense(std))
    ds["runoff_onset_median_n"] = (dims_static, dense(n, np.int32, 0))
    nmad = sums["n_mad"].to_numpy()
    mean, std = _mean_std(sums["s_mad"].to_numpy(), sums["ss_mad"].to_numpy(), np.where(nmad > 0, nmad, np.nan))
    ds["runoff_onset_mad"] = (dims_static, dense(mean))
    ds["runoff_onset_mad_std"] = (dims_static, dense(std))
    ds["runoff_onset_mad_n"] = (dims_static, dense(nmad, np.int32, 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        ntres = sums["n_tres"].to_numpy()
        ds["temporal_resolution_median"] = (dims_static, dense(sums["s_tres"].to_numpy() / np.where(ntres > 0, ntres, np.nan)))
        ds["n_years"] = (dims_static, dense(sums["s_nyears"].to_numpy() / nz))
        ds["latitude_mean"] = (dims_static, dense(sums["s_lat"].to_numpy() / nz))
        nch = sums["n_chili"].to_numpy()
        ds["chili_mean"] = (dims_static, dense(sums["s_chili"].to_numpy() / np.where(nch > 0, nch, np.nan)))
    ds["temporal_resolution_median_n"] = (dims_static, dense(ntres, np.int32, 0))
    ds["terrain_n"] = (dims_static, dense(sums["n_terrain"].to_numpy(), np.int64, 0))

    dims_yearly = dims_static + ("water_year",)
    yr = {k: np.full(shape_static + (len(water_years),), np.nan, dtype=np.float32) for k in ("on", "on_std", "an", "an_std")}
    yr_n = np.zeros(shape_static + (len(water_years),), dtype=np.int32)
    for j, y in enumerate(water_years):
        if f"n_{y}" not in sums:
            continue
        ny = sums[f"n_{y}"].to_numpy(); nyz = np.where(ny > 0, ny, np.nan)
        m, sd = _mean_std(sums[f"s_{y}"].to_numpy(), sums[f"ss_{y}"].to_numpy(), nyz)
        yr["on"][idx + (j,)] = m; yr["on_std"][idx + (j,)] = sd
        m, sd = _mean_std(sums[f"s_an_{y}"].to_numpy(), sums[f"ss_an_{y}"].to_numpy(), nyz)
        yr["an"][idx + (j,)] = m; yr["an_std"][idx + (j,)] = sd
        yr_n[idx + (j,)] = ny
    ds["runoff_onset"] = (dims_yearly, yr["on"])
    ds["runoff_onset_std"] = (dims_yearly, yr["on_std"])
    ds["runoff_onset_n"] = (dims_yearly, yr_n)
    ds["runoff_onset_anomaly"] = (dims_yearly, yr["an"])
    ds["runoff_onset_anomaly_std"] = (dims_yearly, yr["an_std"])
    for x in ("chili", "fcf"):
        r = _pearson(sums[f"n_{x}"].to_numpy(), sums[f"s_{x}"].to_numpy(), sums[f"s_med_{x}"].to_numpy(),
                     sums[f"s_{x}_med"].to_numpy(), sums[f"ss_{x}"].to_numpy(), sums[f"ss_med_{x}"].to_numpy())
        ds[f"{x}_corr"] = (dims_static, dense(r))
        ds[f"{x}_corr_n"] = (dims_static, dense(sums[f"n_{x}"].to_numpy(), np.int32, 0))

    for b in bins:
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
    ds["terrain_n"].attrs["long_name"] = "land pixels (ESA WorldCover not water) of the selection, onset or not"
    ds["latitude_mean"].attrs["long_name"] = "mean latitude of the analysed pixels"
    ds["chili_mean"].attrs["long_name"] = "mean CHILI of the analysed pixels"
    ds["chili_corr"].attrs["long_name"] = "Pearson r of CHILI vs pixel median onset, per cell"
    ds["fcf_corr"].attrs["long_name"] = "Pearson r of forest cover fraction vs pixel median onset, per cell"
    ds.attrs.update({
        "unit_type": unit_type, "unit_id_column": spec["id_col"], "partials_schema": PARTIALS_SCHEMA,
        "selection_fcf_class": str([FCF_CLASSES[c] for c in fcf]) if fcf is not None else "all",
        "selection_slope_class": str([SLOPE_CLASSES[c] for c in slope]) if slope is not None else "all",
        "selection_perennial": str([PERENNIAL_CLASSES[c] for c in perennial]) if perennial is not None else "all",
        "analysed_predicates": str(ANALYSED), "terrain_predicates": str(TERRAIN),
        "statistics": "<var> = cell mean of pixel values; <var>_std = cell std; <var>_n = pixel count "
                      "(runoff_onset_anomaly shares runoff_onset_n). Marginals via aggregate.weighted_mean.",
    })
    if basin_level:
        ds.attrs["hydrobasins_level"] = basin_level
        ds.attrs["hydrobasins_level_stored"] = BASIN_LEVEL_STORED
        ds[spec["unit_dim"]].attrs["long_name"] = f"HydroBASINS level-{basin_level} PFAF_ID"
    summed = [b for b in spec["bins"] if b not in bins] + [k for k in CLASS_KEYS if k not in keep]
    ds.attrs["summed_out"] = ",".join(summed)
    return ds


def cube(sparse, group, filter_tag, water_years, **kw):
    """The named default cube ``group`` (:data:`GROUPS`) for the former filter tag ``filter_tag``
    (:data:`FILTER_VIEWS`): the schema-2 cubes without their ``chili_class`` dimension."""
    g = GROUPS[group]
    ds = view(sparse, g["unit_type"], g["bins"], water_years, fcf=FILTER_VIEWS[filter_tag],
              basin_level=g.get("basin_level"), **kw)
    ds.attrs.update({"group": group, "filter_tag": filter_tag})
    return ds


FACING = ("sunny", "shaded")


def facing_mask(sparse, facing, half_width=45.0, min_slope_class=1):
    """Rows of a sparse table on equator-facing (``'sunny'``) or pole-facing (``'shaded'``) slopes: aspect
    within ``half_width`` deg of the equator direction (south in the northern hemisphere, north in the
    southern; the hemisphere from the row's latitude bin, else from the mean latitude of its analysed
    pixels) or of the pole direction, on slope classes >= ``min_slope_class`` (default >= 10 deg). Flat
    pixels (no aspect) are in neither. This replaces the CHILI warm/cool classes of schema 2."""
    if facing not in FACING:
        raise ValueError(f"facing must be one of {FACING}")
    aspect = sparse["aspect"].to_numpy(dtype=float)
    if "latitude" in sparse and sparse["latitude"].notna().any():
        lat = sparse["latitude"].to_numpy(dtype=float)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            lat = sparse["s_lat"].to_numpy(dtype=float) / sparse["n"].to_numpy(dtype=float)
    equator = np.where(lat >= 0, 180.0, 0.0)
    d = np.abs((aspect - equator + 180.0) % 360.0 - 180.0)         # angular distance from the equator direction
    m = (d <= half_width) if facing == "sunny" else (d >= 180.0 - half_width)
    return m & np.isfinite(aspect) & np.isfinite(lat) & (sparse["slope_class"].to_numpy() >= min_slope_class)


def facing_view(sparse, unit_type, bins, water_years, facing, half_width=45.0, min_slope_class=1, **kw):
    """:func:`view` restricted to sun-facing or shade-facing slopes (:func:`facing_mask`); ``bins`` must
    not include aspect (the facing set is summed)."""
    if "aspect" in bins:
        raise ValueError("facing_view sums aspect; do not keep it in bins")
    sub = sparse[facing_mask(sparse, facing, half_width, min_slope_class)]
    sub.attrs = dict(sparse.attrs)
    return view(sub, unit_type, bins, water_years, **kw)


def add_facing(ds, sparse, unit_type, bins, water_years, half_width=45.0, min_slope_class=1, **kw):
    """Add ``sunny_<var>`` / ``shaded_<var>`` (median onset, MAD, their counts) to a cube ``ds`` with the same
    ``bins`` (e.g. the continents cube): the inputs of the sunny-shaded difference."""
    for facing in FACING:
        f = facing_view(sparse, unit_type, bins, water_years, facing, half_width, min_slope_class, **kw)
        f = f.reindex_like(ds[[d for d in ds.dims if d in f.dims]])
        for v in ("runoff_onset_median", "runoff_onset_median_n", "runoff_onset_mad", "runoff_onset_mad_n"):
            ds[f"{facing}_{v}"] = f[v].fillna(0).astype(f[v].dtype) if v.endswith("_n") else f[v]
    ds.attrs["facing_definition"] = (f"sunny = aspect within {half_width} deg of the equator direction, shaded = within "
                                     f"{half_width} deg of the pole direction, slopes >= {10 * min_slope_class} deg")
    return ds


def reduce_partials(partials, group, filter_tag, water_years):
    """Schema-2-compatible entry point for partials already in memory (tests, the CI fixtures):
    ``cube(sparse_from_frame(partials, unit), group, filter_tag, water_years)``."""
    unit_type = GROUPS[group]["unit_type"]
    return cube(sparse_from_frame(partials, unit_type), group, filter_tag, water_years)


def histogram_view(hist_sparse, unit_type, water_years, perennial=None, by_elevation=False):
    """Onset histograms from a histogram sparse table (``reduce_sparse(..., histograms=True)``): a Dataset
    with ``median_count`` (unit, [elevation], onset_bin) and ``onset_count`` (unit, [elevation], onset_bin,
    water_year); ``onset_bin`` = bin start day of water year (5-day bins). FCF <= 50 % analysed pixels."""
    spec = UNIT_TYPES[unit_type]
    sel = hist_sparse if perennial is None else hist_sparse[hist_sparse["perennial"].isin(list(perennial))]
    keys = ["unit_id"] + (["elevation"] if by_elevation else []) + ["onset_bin"]
    cols = ["h_med"] + [f"h_{y}" for y in water_years if f"h_{y}" in sel]
    s = sel.groupby(keys)[cols].sum().reset_index()
    units = np.array(sorted(s["unit_id"].unique()))
    coords = {spec["unit_dim"]: units, "onset_bin": ONSET_BIN_EDGES[:-1].astype(float)}
    if by_elevation:
        coords["elevation"] = bin_centers("elevation")
    dims = (spec["unit_dim"],) + (("elevation",) if by_elevation else ()) + ("onset_bin",)
    shape = tuple(len(coords[d]) for d in dims)
    idx = [pd.Index(units).get_indexer(s["unit_id"])]
    if by_elevation:
        idx.append(pd.Index(coords["elevation"]).get_indexer(s["elevation"]))
    idx.append(pd.Index(coords["onset_bin"]).get_indexer(s["onset_bin"]))
    idx = tuple(idx)
    med = np.zeros(shape, np.int64); med[idx] = s["h_med"].to_numpy()
    yrs = np.zeros(shape + (len(water_years),), np.int64)
    for j, y in enumerate(water_years):
        if f"h_{y}" in s:
            yrs[idx + (j,)] = s[f"h_{y}"].to_numpy()
    ds = xr.Dataset({"median_count": (dims, med), "onset_count": (dims + ("water_year",), yrs)},
                    coords={**coords, "water_year": np.array(water_years)})
    ds["onset_bin"].attrs["long_name"] = "start day of water year of the 5-day bin"
    return ds


# ---------------------------------------------------------------------------
# read-side helpers (the notebooks' vocabulary for the unified schema)

# quantities whose counts live under another variable's name
_N_ALIAS = {"runoff_onset_anomaly": "runoff_onset_n",
            "runoff_onset_anomaly_std": "runoff_onset_n",
            "runoff_onset_std": "runoff_onset_n",
            "runoff_onset_median_std": "runoff_onset_median_n",
            "runoff_onset_mad_std": "runoff_onset_mad_n",
            "n_years": "runoff_onset_median_n",
            "latitude_mean": "runoff_onset_median_n",
            "chili_mean": "chili_corr_n"}


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
    Schema-3 cubes have no ``chili_class`` dimension, so the default call
    returns the cube unchanged (old notebook code keeps working)."""
    if dim not in ds.dims:
        return ds
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
    """Mirror the partials parquets of ``config.version`` into ``cache_dir``
    (``paths.partials_cache(version)``): download tiles missing from the cache or whose blob changed
    on Azure (a re-mapped tile; compared by etag, recorded in ``.sync_etags.json``, falling back to
    the size for files cached before the manifest existed) and delete cached files Azure no longer
    lists. Re-runs fetch only what changed. Needs the Azure SAS token. Returns the sorted local paths."""
    import json
    import os
    from concurrent.futures import ThreadPoolExecutor
    fs = config.azure_blob_fs
    prefix = f"{settings.PARTIALS_PREFIX}/{config.version}"
    remote = {d['name'].rsplit('/', 1)[-1]: d for d in fs.ls(prefix, detail=True) if d['name'].endswith('.parquet')}
    manifest_path = cache_dir / '.sync_etags.json'
    etags = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    local = {f for f in os.listdir(cache_dir) if f.endswith('.parquet')}
    stale = local - set(remote)
    for f in stale:
        os.remove(cache_dir / f)
        etags.pop(f, None)

    def changed(name):
        if name not in local:
            return True
        if name in etags:
            return etags[name] != remote[name].get('etag')
        return os.path.getsize(cache_dir / name) != remote[name].get('size')

    todo = sorted(n for n in remote if changed(n))
    n_new = sum(n not in local for n in todo)
    log(f"partials on Azure: {len(remote)} | cached: {len(local) - len(stale)} | "
        f"downloading {len(todo)} ({n_new} new, {len(todo) - n_new} changed) | dropped stale {len(stale)}")

    def fetch(name):
        fs.get_file(remote[name]['name'], str(cache_dir / f"{name}.part"))
        os.replace(cache_dir / f"{name}.part", cache_dir / name)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(fetch, todo))
    etags.update({n: remote[n].get('etag') for n in remote})
    manifest_path.write_text(json.dumps(etags))
    return sorted(cache_dir / f for f in remote)
