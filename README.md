# global_snowmelt_runoff_onset_analysis

Analyses of the [global snowmelt runoff onset dataset](https://doi.org/10.5281/zenodo.16953614)
([Gagliano et al., 2026, *Earth System Science Data*](https://doi.org/10.5194/essd-18-5871-2026)):
a map/reduce pipeline that aggregates the 80 m dataset into per-unit statistics (mountain ranges,
river basins, continents) and the notebooks that analyse them — elevation and aspect structure of
runoff onset, its interannual anomalies, its sensitivity to spring temperature, basin-scale summaries,
a single-range worked example.

Dataset creation, evaluation and the interactive map live in the companion repo
[egagli/global_snowmelt_runoff_onset](https://github.com/egagli/global_snowmelt_runoff_onset), whose
package this one imports for the dataset store, the configuration and the shared plot styling.

## Setup

The production repo is cloned **side-by-side** (the pixi environment editable-installs it from
`../global_snowmelt_runoff_onset`; `Config` resolves the config file, credentials and tile registries
inside that clone). The workflows and the CI check it out at the commit named in `pixi.toml`.

```bash
git clone https://github.com/egagli/global_snowmelt_runoff_onset.git
git clone https://github.com/egagli/global_snowmelt_runoff_onset_analysis.git
cd global_snowmelt_runoff_onset_analysis
pixi install                # the default environment (pipeline + analysis + jupyter)
pixi run setup-kernel       # registers the 'gsro-analysis (pixi)' jupyter kernel
pixi run lab
```

Credentials, all resolved automatically:

| Credential | Needed for | Source |
| --- | --- | --- |
| Azure SAS token | the dataset store and every Azure product (the ancillary and ERA5-Land repos, partials, mirrored cubes) | `../global_snowmelt_runoff_onset/config/sas_token.txt` or `AZURE_STORAGE_SAS_TOKEN` |
| Earth Engine service key | CHILI (ancillary build), ERA5-Land acquisition, the GTOPO30 histogram rebuild, `river_basins/snow_water` | `../global_snowmelt_runoff_onset/config/ee_key.json`, via `settings.initialize_earthengine()` (works headless) |
| none | the public multiscale pyramid, every notebook that reads only the local cubes and tables, the Evan & Eisenman comparison | — |

Compute is local or GitHub Actions; there is no cluster service.

## Layout

```text
gsro_analysis/        the package (table below)
pipeline/             scripts/ (the two wave workers, the dispatcher, the ERA5-Land store, zonal join, reduce, metrics),
                      pipeline.ipynb (state, reset, work list + dispatch, one tile end to end, local
                      stages), tile_data/ (the fleet work list); README.md = mechanics
analyses/<unit>/      one folder per aggregation unit (global = continents cube, mountain_ranges,
                      river_basins) + case_studies/sierra_nevada and climate/; figures/<version>/,
                      results/<version>/ tracked per dataset version; README.md = notebook table
data/                 shared inputs (README.md says where each comes from) — gitignored except the
                      40 KB GTOPO30 histogram
aggregated_results/   pipeline outputs, version-first — gitignored, mirrored to Azure
  <version>/<group>/all_<group>_<filter>.nc    the cubes (mountain_ranges, river_basins, continents;
                                               river_basins_l6 and continents_aspect on demand)
  <version>/era5_zonal/                        per-unit ERA5-Land anomaly zonal means
  <version>/partials/                          local cache of the per-tile partial sums
tests/                fixture partials for the credential-free CI smoke test
.github/workflows/    get_era5_land_data.yml, get_ancillary_data.yml, process_tiles.yml (workflow_dispatch), ci.yml
```

**Version discipline.** Every output path names the dataset version it came from (`config.version`,
e.g. `v10`). The version is chosen in one place, `gsro_analysis.settings.CONFIG_FILE` (override per
shell with `GSRO_CONFIG=…`), and every notebook and script builds its `Config` through
`settings.load_config()`. Pixel filters are named predicate sets (`aggregate.FILTERS`: `fcf_lte_50`,
the rule of every analysis, and `full_dataset`), applied in the fleet's map step and recorded in every
file name and attrs. Water years come from `config.water_years`, never literals. `GSRO_AGGREGATED_ROOT`
and `GSRO_OUTPUT_ROOT` redirect the cubes and the figures/results to a test set.

## Pipeline

Three waves of GitHub Actions workflows, dispatched by hand in this order and each re-dispatched
until its plan job reports nothing remaining, then the local stages:

```text
 1  Get ERA5-Land data      (get_era5_land_data.yml -> scripts/era5_land.py; one job per water year)
      ERA5-Land monthly (Earth Engine via xee, native 0.1 deg grid) -> ONE icechunk repo per version
      snowmelt_runoff_onset_analysis/era5_land/<version>/era5_land: one commit per water year + the anomaly group

 2  Get ancillary data      (get_ancillary_data.yml -> scripts/build_ancillary_batch.py; ~2.5 min per tile)
      Cop-DEM, CHILI (Earth Engine), snow class, WorldCover, forest cover, GMBA / HydroBASINS level-6 /
      continent ids -> the tile's window of the DATASET grid (EPSG:4326, 0.00072 deg) -> ONE icechunk repo per grid
      snowmelt_runoff_onset_analysis/ancillary/<version>_grid/ancillary: compact ints, one chunk and one commit per tile

 3  Process tiles to parquets (process_tiles.yml -> scripts/process_tiles_batch.py; ~1 min per tile, no Earth Engine)
      store window + ancillary window (same pixel indices) -> the tile's 80 m UTM grid (onset/DEM/CHILI/forest cover
      bilinear, categorical and id layers nearest) -> slope + aspect from the UTM DEM -> pixel table (in memory)
      -> PARTIAL SUMS per (filter, unit, bins, CHILI class) -> snowmelt_runoff_onset_analysis/partials/<version>/tile_RRR_CCC.parquet
      [keep_pixels] -> the per-pixel parquet (opt-in)

 local (minutes each; pixi run era5-zonal / reduce / metrics)
   scripts/era5_zonal.py       anomaly group x (GMBA, HydroBASINS) polygons -> era5_zonal/*.nc
   scripts/reduce_partials.py  partials -> the cubes (+ ERA5 zonal means, GTOPO30, names)
   scripts/range_metrics.py    mountain-range cube -> results/<version>/mountain_range_metrics.csv
   analyses/<unit>/*.ipynb     notebooks -> figures/<version>/
```

The two stores are ledgers: a tile or a water year is done when it has a commit, a failed job commits
nothing, and each plan job folds the commit history to list what remains (the production repo's
icechunk fleet pattern). Every statistic the analyses read is reducible (means, standard deviations and
correlations from sums), so a tile job emits sums instead of a pixel table and the reduce is a laptop
job. The map keys every unit type by the finest axes it will ever need (basins at HydroBASINS level 6,
continents by latitude x elevation x aspect); the reduce derives the coarser default cubes (level-5
basins, latitude x elevation continents) exactly, and the finer ones on demand. Tabulating on the UTM
grid keeps every row a ~6,400 m² pixel, so pixel counts are area without weights. Mechanics, ledgers,
what a partials row is, sizing, redo recipes: [pipeline/README.md](pipeline/README.md). Lineage of
every input from source to analysis: [docs/aggregation_lineage.md](docs/aggregation_lineage.md).

### The aggregate schema

One schema for every unit type (`aggregate.UNIT_TYPES` = the map keys, `aggregate.GROUPS` = the cubes):

| cube | dims | extras merged by the reduce |
| --- | --- | --- |
| `mountain_ranges` | `mountain_range × elevation(100 m) × aspect(15°) × chili_class × water_year` | GMBA_V2_ID, centroid, continent; ERA5-Land anomaly zonal means `(range, water_year, month)` |
| `river_basins` (HydroBASINS level 5, derived from the stored level-6 ids) | `river_basin × elevation × chili_class × water_year` | ERA5-Land zonal means `(basin, water_year, month)` |
| `continents` | `continent × latitude(1°) × elevation × chili_class × water_year` | `dem_pixel_count` (GTOPO30 land histogram) |
| on demand: `river_basins_l6`, `continents_aspect` | level-6 basins; `continent × latitude × elevation × aspect × …` | as above |

Variables: `runoff_onset_median`, `runoff_onset_mad`, `temporal_resolution_median` (static),
`runoff_onset`, `runoff_onset_anomaly` (per water year), each as the bin **mean** with `<var>_std` and
`<var>_n`; `n_years`; `chili_corr`, `fcf_corr`. Nothing is precomputed: whole-basin means, elevation
profiles and range mean anomalies are `aggregate.weighted_mean(ds, var, dims)`; `aggregate.collapse`
folds the CHILI axis exactly; `aggregate.threshold` masks thin bins; `stats.prepare_mountain_ranges`
applies the analyses' rules in one call.

## Running a dataset version

1. **1. Get ERA5-Land data** (the Actions tab lists the three workflows in this order): no inputs needed; it creates the version's store if missing, fetches only the
   water years without a commit, then builds the anomaly group (`start_fresh` deletes the store first).
2. **2. Get ancillary data**: creates the grid generation's ancillary store if missing and builds every tile
   without a commit (`start_fresh` deletes the store first). Re-dispatch until the plan job reports 0.
3. **3. Process tiles to parquets**: maps every tile that has an ancillary commit and no partials blob
   (`start_fresh` deletes the partials and pixel tables first). Re-dispatch until 0 remaining. Each job
   logs one line per step per tile.
4. **Local stages**: `pixi run era5-zonal`, `pixi run reduce` (add `-- --mirror` to copy the cubes to
   Azure), `pixi run metrics`.
5. **Notebooks**: see [analyses/README.md](analyses/README.md); the two composite world maps build their
   own inset sweeps.

`pipeline/pipeline.ipynb` has a cell for each step (state, reset, work list, dispatch, one tile end to
end, local stages). To redo a tile, delete its partials blob and re-dispatch wave 3; recommitting a tile
in wave 2 supersedes its ancillary. A change of a unit definition alone is
`datacube.refresh_unit_layers` (one commit per tile, no Earth Engine) plus a re-map.

## The `gsro_analysis` package

| Module | Contents |
| --- | --- |
| [paths.py](gsro_analysis/paths.py) | version-aware, repo-root-anchored paths: `figdir`, `resultsdir`, `aggregate`, `era5_zonal`, `partials_cache`, `gtopo30_histogram`; the two environment overrides |
| [settings.py](gsro_analysis/settings.py) | the one place the dataset version is chosen; Azure prefixes; external source URLs (GMBA, BasinATLAS, continents, forest cover, snow class); Earth Engine init; `cached_source` |
| [ledger.py](gsro_analysis/ledger.py) | the icechunk ledger shared by the two stores: storage handles, `commit_with_retry`, the ancestry fold `commit_records` |
| [datacube.py](gsro_analysis/datacube.py) | the ancillary store on the dataset grid (`initialize_ancillary_store`, `build_ancillary_window`, `write_ancillary_tile`, `open_ancillary_window`, `refresh_unit_layers`), the UTM process step (`tabulate_tile`, `process_tile`), partials/pixel-table writers, `reset_version` (dry run by default) |
| [aggregate.py](gsro_analysis/aggregate.py) | `FILTERS`, `UNIT_TYPES`, `GROUPS`; the map `tile_partials` and the reduce `reduce_partials`; `open_aggregate`, `weighted_mean`, `collapse`, `threshold`, `elevation_relative`; range names and metadata; the GTOPO30 histogram |
| [era5.py](gsro_analysis/era5.py) | the ERA5-Land icechunk store: template, native-grid acquisition (`acquire_water_year`), commit ledger (`status`), `build_anomaly`, `open_era5_land`, `open_anomaly`, the zonal join `zonal_anomalies`, Equal Earth on the fly, `pixelwise_correlations` |
| [stats.py](gsro_analysis/stats.py) | `prepare_mountain_ranges`, `basin_summary`, two lapse-rate definitions, `spring_temperature_sensitivity`, `climate_regressions`, `range_metrics_gdf` |
| [results.py](gsro_analysis/results.py) | `save_result_table`: the production helper plus `_analysis_git_sha` |
| [plotting.py](gsro_analysis/plotting.py), [colorbars.py](gsro_analysis/colorbars.py) | polar triplets and anomaly panels; annotated colorbars as presets over the production `plot_utils` |
| [world_maps.py](gsro_analysis/world_maps.py) | the two inset world maps (page model in mm, label blocks from `analyses/mountain_ranges/label_layout.csv`, choropleth ramps, `layout_report`) |

## Conventions

Statistics come from the cube via the read-side helpers, never re-derived by hand; every quoted number
exists in a `results/<version>/` table stamped with `_git_sha` (production package), `_analysis_git_sha`
(this repo) and `_written_at`; figures are PNG at ≤ 300 dpi, one tracked set per dataset version, sweeps
gitignored; figure styling follows the production `plot_utils` conventions; notebook outputs are
stripped by `nbstripout`. Nothing on Azure is deleted by the code except `datacube.reset_version`, which
is a dry run until confirmed.

License: MIT (`LICENSE`). When using these analyses, cite the dataset publication (`CITATION.cff`):

> Gagliano, E., Shean, D., and Henderson, S.: A global high-resolution dataset of snowmelt runoff onset
> timing from Sentinel-1 SAR, 2015–2024, Earth Syst. Sci. Data, 18, 5871–5894,
> <https://doi.org/10.5194/essd-18-5871-2026>, 2026.
