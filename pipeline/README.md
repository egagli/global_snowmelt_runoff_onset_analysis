# pipeline

Three waves of GitHub Actions workflows over the tiles, then local stages. The engine lives in
[`gsro_analysis/datacube.py`](../gsro_analysis/datacube.py) (the ancillary store and the process step),
[`gsro_analysis/era5.py`](../gsro_analysis/era5.py) (the ERA5-Land store),
[`gsro_analysis/ledger.py`](../gsro_analysis/ledger.py) (the shared icechunk ledger) and
[`gsro_analysis/aggregate.py`](../gsro_analysis/aggregate.py) (the map `tile_partials`, the reduce
`reduce_partials`); `scripts/` are the entry points; `pipeline.ipynb` is the operator's notebook (state,
reset, work list and dispatch, one tile end to end, local stages). Lineage of every input from source to
analysis, and the design notes: [docs/aggregation_lineage.md](../docs/aggregation_lineage.md).

| Wave / stage | Where | Entry point | Reads | Writes |
| --- | --- | --- | --- | --- |
| 1 Get ERA5-Land data (per **version**) | workflow, one job per water year | `scripts/era5_land.py` | ERA5-Land monthly (Earth Engine via xee, native 0.1° grid) | `snowmelt/snowmelt_runoff_onset_analysis/era5_land/<version>/era5_land`: ONE icechunk repo, the 8 variables on `(water_year, month, latitude, longitude)`, one commit per water year, plus the `anomaly` group (one commit) |
| 2 Get ancillary data (once per **grid**) | workflow, ~2.5 min per tile | `scripts/build_ancillary_batch.py` → `datacube.build_ancillary_window` | Cop-DEM, CHILI (EE), snow class, WorldCover, forest cover, GMBA, BasinATLAS level 6, continents | `snowmelt/snowmelt_runoff_onset_analysis/ancillary/<version>_grid/ancillary`: ONE icechunk repo on the dataset grid (EPSG:4326, 0.00072°), 8 layers with compact int encodings, one chunk and one commit per tile |
| 3 Process tiles to parquets (per **version**) | workflow, ~1 min per tile, no EE | `scripts/process_tiles_batch.py` → `datacube.process_tile` | the dataset store's tile window + the ancillary store's tile window | both reprojected onto the tile's 80 m UTM grid, slope + aspect derived, pixel table in memory → `snowmelt/snowmelt_runoff_onset_analysis/partials/<version>/tile_RRR_CCC.parquet` (partial sums for both filter tags × three unit types, ~0.1–1 MB); `keep_pixels`: also the pixel table `…/parquets/<version>/tile_RRR_CCC.parquet` (~35 MB) |
| 4 ERA5 zonal | local | `scripts/era5_zonal.py` | anomaly group, GMBA + BasinATLAS polygons (level 5; `--units river_basins_l6` for level 6), pyramid masks | `aggregated_results/<version>/era5_zonal/era5_anomaly_<unit_type>.nc` |
| 4 reduce | local | `scripts/reduce_partials.py` | the partials (auto-cached locally), GMBA/continents, `data/gtopo30_lat_elev_histogram.nc`, the zonal files | `aggregated_results/<version>/<group>/all_<group>_<filter>.nc` (`--groups` adds `river_basins_l6`, `continents_aspect`; `--mirror` → `snowmelt/snowmelt_runoff_onset_analysis/aggregated/<version>/`) |
| 5 metrics | local | `scripts/range_metrics.py` | the mountain-range cube | `analyses/mountain_ranges/results/<version>/mountain_range_metrics.csv` |

## What a "partials" row is

A partial sum is one tile's contribution to one cell of the final cube. The map
(`aggregate.tile_partials`) groups the tile's pixels by (filter, unit type, unit id, elevation bin,
aspect bin, latitude bin, CHILI class) — the axes each unit type uses (`aggregate.UNIT_TYPES`) — and
for each group writes the pixel count and the sums the statistics need — Σ median, Σ median², Σ onset
per year, Σ anomaly per year, Σ CHILI·median, … (83 columns). A range that spans 12 tiles therefore
shows up as 12 rows per (bin, class); the reduce (`aggregate.reduce_partials`) adds the rows and turns
them into means (Σx/n), standard deviations (√(Σx²/n − mean²)) and correlations. That is why the
partials are "partial": each tile only knows its share of a range, basin or continent, and only their
sum is the answer.

## Why partial sums

Every statistic the analyses read is reducible — bin means and standard deviations from
(n, Σx, Σx²), the CHILI/FCF correlations from (n, Σx, Σy, Σxy, Σx², Σy²) — and no notebook ever read
a bin-level median. So the fleet writes sums, not pixels: the tile products shrink from ~50 MB to
<1 MB, stage 2 becomes a pandas reduce that runs in minutes on a laptop, and every unit type is
produced in one pass. Validated 2026-08-26 and again 2026-09-03 (after the schema change below) on the
five v10 dry-run tiles: the reduce reproduces an independent pandas computation to float32 precision
(counts exactly) for every group and filter.

## Units, levels and groups (decided 2026-09-02, in the code since 2026-09-03)

The map keys every unit type by the **finest axes it will ever need**, because a coarser axis is a
free `groupby` at reduce time while a finer one is a fleet re-map:

- **River basins** are stored at HydroBASINS **level 6** (`settings.BASIN_ATLAS_LAYER`,
  six-digit `PFAF_ID`). Pfafstetter codes nest by digit prefix, so the default `river_basins` cube
  (level 5, `PFAF_ID // 10`) is exact; `--groups river_basins_l6` writes the level-6 cube (run
  `era5_zonal.py --units river_basins_l6` first for its ERA5 merge). Any level ≤ 6 is one more
  entry in `aggregate.GROUPS`.
- **Continents** carry the 24-bin **aspect** key. The default `continents` cube is latitude ×
  elevation (aspect summed out, so flat pixels with no aspect still count); `--groups
  continents_aspect` writes the latitude × elevation × aspect cube (drops the ~0.06 % flat pixels).
- **CHILI** is the fourth key of every unit type; `aggregate.collapse` folds it exactly on read.
- **Filters** (`aggregate.FILTERS`): `fcf_lte_50` = seasonal snow class ≠ 4 (ephemeral) + WorldCover
  ∉ {50, 80} + 0 ≤ fcf ≤ 50 (the ≥ 0 bound also drops nodata); `full_dataset` = the first two only.

Changing a unit definition after the ancillary exists is `datacube.refresh_unit_layers(config, row,
col)` (re-rasterizes the id layers of a stored tile into the ancillary store, one commit, no Earth
Engine) followed by a re-map of the tile (`process_tile`).

## Wave mechanics

The ~4,320-tile campaign runs as GitHub Actions matrices (`.github/workflows/get_ancillary_data.yml`,
`process_tiles.yml`), the production repo's icechunk fleet pattern:

- **ledgers**: wave 2 = a commit with metadata (`kind: ancillary_tile, tile: [row, col]`) in the
  ancillary repository, folded from its history (`datacube.completed_ancillary_tiles`); wave 3 = the
  partials parquet exists (single blob = atomic, one list call). A failed tile commits or writes nothing
  and is re-listed by the next dispatch.
- **dispatch**: `scripts/get_remaining_work.py --stage ancillary|partials` folds the ledgers, diffs
  against `tile_data/ancillary_tiles_v10.txt` (every tile whose composites hold data, from the dataset's
  commit history; `pipeline.ipynb` regenerates it) and writes a batch manifest that every matrix job of
  the run reads; wave 3 dispatches only tiles WITH an ancillary commit. The workers
  (`build_ancillary_batch.py`, `process_tiles_batch.py`) isolate errors per tile and log one line per
  step (target grid, every layer with its source and timing, the write and commit; the two window
  reads, the UTM reprojection and tabulation, the map). Re-dispatch until the plan job reports 0.
- **start_fresh** (off by default) deletes that wave's products of the version first
  (`datacube.reset_version`: the ancillary repository, or the partials and pixel tables) — the only
  deletion the workflows can make.
- **hosting**: the GitHub-hosted runners of this repository with the repository secrets
  (`AZURE_STORAGE_SAS_TOKEN`, `AZURE_STORAGE_ACCOUNT`, `EE_SERVICE_ACCOUNT_KEY`; wave 3 needs no
  Earth Engine); the production repo is checked out at the commit pinned in the workflow `env`
  (mirrored in `pixi.toml`). One run per workflow and config file at a time (concurrency groups).
- **sizing**: wave 2 ~2.5 min per tile (the forest-cover read dominates) plus ~12 min fixed per job
  (env + the 2.7 GB BasinATLAS cache warm) → 36 tiles/job; wave 3 ~1 min per tile, ~3 min fixed →
  72 tiles/job. Failure policy in `build_ancillary_window`: a raster layer fills with nodata only
  outside its documented latitude coverage (`datacube.LAYER_LAT_COVERAGE`) or on the source's own
  no-data signal; any other failure raises (failure = no output, never wrong output).
- **grids**: the ancillary lives on the dataset grid (one grid generation is built with one
  `pixi.lock`: rebuilding with different GDAL/PROJ versions changes the terrain layers at the metre
  level); the tabulation grid is the tile's UTM 80 m grid, onto which both windows are reprojected
  (`datacube.UTM_RESAMPLING`), so pixel counts stay area. The design note in
  [docs/aggregation_lineage.md](../docs/aggregation_lineage.md) has the trade-offs.
- **memory**: wave 2 peaks at ~3 GB per tile, wave 3 at ~3–4 GB (the tabulate); fine on the 16 GB
  public runners.

## The ERA5-Land store

One icechunk repository per dataset version (`gsro_analysis/era5.py`), the production repo's fleet
pattern applied to eleven work units: the `Get ERA5-Land data` workflow's plan job creates the repository
with an empty, all-NaN template if it does not exist (every water year of the config, 9 hemisphere-aware
months, the native 1800 × 3600 grid with cell centres at multiples of 0.1°, north-down), folds the
commit history to list the water years without a commit, and runs one job per missing year; each job
fetches its year from Earth Engine (`ee.data.computePixels`, one float32 band per half-globe request, no
resampling), writes its chunks (one chunk = one variable × year × month, so concurrent jobs never touch
the same chunk) and makes ONE commit with QA metadata. The anomaly job then builds the `anomaly` group
(every variable minus its per-pixel median over all water years) once every year is committed, and
skips when it is current; a year committed after the anomaly marks it stale. A failed job commits
nothing; re-dispatch until nothing remains. Acquisitions are never copied between versions; a new
version re-acquires everything (~5 min per year in parallel, ~20 min for the anomaly). Then
`scripts/era5_zonal.py --overwrite` and `scripts/reduce_partials.py` rerun. `start_fresh` (off by
default) deletes the version's repository first — the only deletion the workflow can make.

## Validation aids

The v9 per-tile parquets (`settings.V9_TILE_PARQUET_PREFIX`) are the reference for the tile-for-tile
comparison in `pipeline.ipynb` (v9 tile (r, c) == v10 tile (r+2, c); v9 basin ids are level 5, i.e.
`PFAF_ID // 10` of the current ones). `datacube.partials_from_pixel_table` backfills partials for a
tile whose pixel table exists. `tests/` holds two fixture partials for the credential-free CI reduce.
