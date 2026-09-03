# pipeline

Map/reduce over tiles. The engine lives in [`gsro_analysis/datacube.py`](../gsro_analysis/datacube.py)
(stage 0/1: ancillary, tabulate, `process_tile`) and [`gsro_analysis/aggregate.py`](../gsro_analysis/aggregate.py)
(the map `tile_partials`, the reduce `reduce_partials`); `scripts/` are the batch entry points;
`pipeline.ipynb` is the operator's notebook (state, reset, work list and dispatch, one tile end to end,
local stages). Lineage of every input from source to analysis, and the design notes behind the unit
levels and the CHILI axis: [docs/aggregation_lineage.md](../docs/aggregation_lineage.md).

| Stage | Where | Entry point | Reads | Writes |
| --- | --- | --- | --- | --- |
| 0 ancillary (once per **grid**) | fleet job | `scripts/run_fleet_batch.py` | Cop-DEM, CHILI (EE), snow class, WorldCover, FCF, GMBA, BasinATLAS level 6, continents | `snowmelt/analysis/ancillary/<version>_grid/tile_RRR_CCC.zarr` (compact int encodings, ~8 MB) + `_complete/tile_RRR_CCC.json` |
| 1 tabulate + map (per dataset **version**) | same fleet job | `scripts/run_fleet_batch.py` → `datacube.process_tile` | icechunk store tile window + the ancillary zarr | `snowmelt/analysis/partials/<version>/tile_RRR_CCC.parquet` — partial sums for both filter tags × three unit types (~0.1–1 MB). `--keep-pixels`: also the pixel table `snowmelt/analysis/parquets/<version>/tile_RRR_CCC.parquet` (~35 MB) |
| ERA5 | `ERA5 Acquire` workflow (16 GB runner) | `scripts/fetch_era5.py` | ERA5-Land monthly (EE) | `snowmelt/analysis/era5_data/<version>/…` (year stores + anomaly store, verified + marked) |
| 2 ERA5 zonal | local | `scripts/era5_zonal.py` | anomaly store, GMBA + BasinATLAS polygons (level 5; `--units river_basins_l6` for level 6), pyramid masks | `aggregated_results/<version>/era5_zonal/era5_anomaly_<unit_type>.nc` |
| 2 reduce | local | `scripts/reduce_partials.py` | the partials (auto-cached locally), GMBA/continents, `data/gtopo30_lat_elev_histogram.nc`, the zonal files | `aggregated_results/<version>/<group>/all_<group>_<filter>.nc` (`--groups` adds `river_basins_l6`, `continents_aspect`; `--mirror` → `snowmelt/analysis/aggregated/<version>/`) |
| 3 metrics | local | `scripts/range_metrics.py` | the mountain-range cube | `analyses/mountain_ranges/results/<version>/mountain_range_metrics.csv` |

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
col)` (stage 0b: re-rasterizes the id layers of a stored tile in place, no Earth Engine, nothing
deleted) followed by a re-map of the tile (`process_tile`); the five dry-run tiles went through
exactly that on 2026-09-03.

## Fleet mechanics

The full ~4,320-tile campaign runs as a GitHub Actions fleet (`.github/workflows/pipeline_fleet.yml`),
the production repo's icechunk fleet pattern adapted to blob-existence ledgers:

- **done-check**: stage 0 = marker blob in `ancillary/<grid>/_complete/` (flat dir, written only
  AFTER a verified re-read — zarr metadata appears at the START of a write, so its presence
  can't be trusted). Stage 1 = the partials parquet exists (single blob = atomic).
- **dispatch**: `scripts/get_remaining_work.py` lists both ledgers (2 list calls), diffs against
  `tile_data/ancillary_tiles_v10.txt` (every tile whose composites hold data, from the icechunk
  history; `pipeline.ipynb` regenerates it), writes a batch manifest; `scripts/run_fleet_batch.py`
  runs `process_tile` per tile with per-tile error isolation (a failed tile leaves NO output and
  is re-listed next dispatch). Re-dispatch until the plan job reports 0 remaining.
- **hosting**: the workflows run on the GitHub-hosted runners of this repository with the three
  repository secrets (`AZURE_STORAGE_SAS_TOKEN`, `AZURE_STORAGE_ACCOUNT`, `EE_SERVICE_ACCOUNT_KEY`);
  the production repo is checked out at the commit pinned in the workflow `env` (mirrored in
  `pixi.toml`). A concurrency group keeps one fleet run per config file.
- **sizing** (2026-08-25 canary): ~2.5 min per tile plus ~12 min fixed per job (env + the 2.7 GB
  BasinATLAS cache warm) → 36 tiles/job (~1.7 h, ~120 jobs, ~200 job-hours). Failure policy in
  `build_ancillary_tile`: a raster layer fills with nodata only outside its documented latitude
  coverage (`datacube.LAYER_LAT_COVERAGE`) or on the source's own no-data signal; any other
  failure raises (failure = no output, never wrong output).
- **memory**: the tabulate step peaks at ~3–4 GB per tile; fine on the 16 GB public runners.

## ERA5 for a new version

The per-water-year stores are version-independent acquisitions: either server-side copy the
existing years to `era5_data/<version>/` and fetch only the new year(s) (`ERA5 Acquire` with
`copy_from_version`), or refetch every year from Earth Engine; then build the anomaly store (base
period = all water years in the stack, recorded in attrs), then rerun `scripts/era5_zonal.py` and
`scripts/reduce_partials.py`.

## Validation aids

The v9 per-tile parquets (`settings.V9_TILE_PARQUET_PREFIX`) are the reference for the tile-for-tile
comparison in `pipeline.ipynb` (v9 tile (r, c) == v10 tile (r+2, c); v9 basin ids are level 5, i.e.
`PFAF_ID // 10` of the current ones). `datacube.partials_from_pixel_table` backfills partials for a
tile whose pixel table exists. `tests/` holds two fixture partials for the credential-free CI reduce.
