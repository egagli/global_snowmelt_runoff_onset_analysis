# analyses

One folder per aggregation unit (`continents/`, `mountain_ranges/`, `river_basins/`), plus the Sierra Nevada case
study (reads the dataset store directly) and `climate/` (ERA5-Land products and a comparison with a published
sensitivity model). Every unit folder has the same shape:

```text
analyses/<unit>/
  0_aggregate_by_<unit>.ipynb      the fleet's partial sums -> the unit's cube(s), ERA5-Land zonal means, the unit's metrics table
  <topic>.ipynb                    the analyses: they read the cube and the metrics table, and write figures
  data/geometries/                 the unit's polygons, downloaded once (gitignored; the small tracked tables are listed below)
  data/aggregation/<version>/      all_<unit>_<filter>.nc (one cube per pixel filter), era5_anomaly_<unit>.nc (gitignored, regenerable)
  results/<version>/               the tables (tracked): <unit>_metrics.csv and what the topic notebooks write
  figures/<version>/               the figures (tracked, except the machine-generated sweeps)
```

Run the `0_aggregate` notebook of a unit first (after the three GitHub Actions workflows), then any topic notebook: each
opens with a markdown cell that says what it reads and writes; run its setup cells (imports, `config =
settings.load_config()`, the reads) first, then any section. Display-range lists are filtered to what exists in the cube,
so a partially processed dataset version runs.

## The notebooks

| Folder | Notebook | Reads | Writes (`data/aggregation/<version>/`, `results/<version>/`, `figures/<version>/`) |
| --- | --- | --- | --- |
| `continents/` | `0_aggregate_by_continent` | `partials/<version>/` (downloaded from Azure with the SAS token, else the local cache), `data/gtopo30_lat_elev_histogram.nc` (tracked) | `all_continents_<filter>.nc`; `all_continents_aspect_<filter>.nc` on request |
| `continents/` | `lat_elev_binning` | the continents cube (`fcf_lte_50`; `full_dataset` for the share of pixels above 5000 m) | `continent_lat_and_elev_binned.png` (median, MAD, sunny–shaded difference per continent, annotated colorbars drawn in place), `continent_lat_and_elev_binned_anomaly.png`, FCF/CHILI correlation panels |
| `continents/` | `sunny_shaded` | the continents cube | `global_scatter_latitude_vs_onset.png` (seasonal modulation of the sunny–shaded offset), offset-vs-onset and MAD-vs-onset scatters |
| `continents/` | `lapse_rates` | the continents cube | `lapse_rates_by_latitude_groupby_continents.png`, `global_lapse_rates_by_continent.png`; latitude + elevation regressions |
| `mountain_ranges/` | `0_aggregate_by_mountain_range` | `partials/<version>/`, GMBA v2 + USGS continents (web), the ERA5-Land anomaly group (Azure) | `all_mountain_ranges_<filter>.nc`, `era5_anomaly_mountain_ranges.nc`, `mountain_range_metrics.csv` |
| `mountain_ranges/` | `topography` | the mountain-range cube, the metrics table | `polar_triplets_all_ranges.png`, `global_lapse_rates_by_mountain_range.png`, `global_mad_by_mountain_range.png` |
| `mountain_ranges/` | `topography_triplet_composite_figure` | the mountain-range cube, the metrics table joined to the GMBA polygons, `label_layout.csv`, the hillshade | the per-range triplet sweep `triplets/` (gitignored), `polar_triplet_legend.png`, the composite `global_lapse_rates_with_triplets_map.png` (`gsro_analysis.world_maps`) |
| `mountain_ranges/` | `anomalies` | the mountain-range cube, GMBA polygons | `snowmelt_onset_anomalies_by_mountain_range.png`, `regional_panels/*.png` (gitignored sweep), `regional_panels/combined_annual/*.png`, `median_runoff_and_anomaly_subset.png` |
| `mountain_ranges/` | `temperature_sensitivity` | the mountain-range cube (with the ERA5-Land zonal means), the metrics table | `snowmelt_onset_spring_temp_sensitivity_histograms.png`, `mountain_range_spring_temp_vs_runoff_onset_anomaly_scatterplots_theil_sen.png`, regression heatmaps |
| `mountain_ranges/` | `spring_temperature_sensitivity_composite_figure` | as the triplet composite, plus the ERA5-Land zonal means | the per-range scatterplot sweep `anomaly_scatterplots/` (gitignored), the composite `global_temperature_sensitivity_map.png` |
| `mountain_ranges/` | `teleconnections` | NOAA PSL Niño 3.4 / PDO / PNA monthly indices (read live) | `teleconnection_indices.png` |
| `river_basins/` | `0_aggregate_by_river_basin` | `partials/<version>/`, BasinATLAS level-5 polygons (cached gdb), `data/geometries/hydrobasins_level6_population.csv` (tracked), the ERA5-Land anomaly group (Azure) | `all_river_basins_<filter>.nc` (HydroBASINS level 5), `era5_anomaly_river_basins.nc`, `river_basin_metrics.csv`; the level-6 cube on request |
| `river_basins/` | `basin_onset` | the river-basin cube, the metrics table, the level-5 polygons, the hillshade, `river_basin_snow_water.csv` if present, the World Bank major rivers (web) | `global_median_runoff_onset_and_mad.png`, `global_runoff_onset_anomaly_by_basin.png` and, for `hma` and `wus`: `<region>_runoff_onset_anomaly_by_basin.png`, `<region>_basin_population_and_pct_precip_as_snow.png`, `<region>_basin_runoff_onset_population_pct_precip_as_snow.png`, `<region>_basin_runoff_onset_vs_pct_precip_as_snow.png`, `<region>_basin_runoff_onset_vs_pct_precip_as_snow_largest_anom.png` |
| `river_basins/` | `snow_water` | the river-basin cube, the metrics table, the level-5 polygons, ERA5-Land via Earth Engine (April-1 / October-1 SWE, precipitation-as-snow share) | `river_basin_snow_water.csv` (read by `basin_onset`), `pct_precip_as_snow.png`, `mad_median_swe.png`, `swe_median_mad_runoff_onset_pop.png`, SWE / water-volume vs onset scatters |
| `case_studies/sierra_nevada/` | `sierra_nevada` | the dataset store (4x coarsened), the ERA5-Land store, the hillshade, forest cover (FCF ≤ 50 % masked like the pipeline); peaks at ~12 GB, so run it alone on a 16 GB machine with `MALLOC_ARENA_MAX=2 MALLOC_MMAP_THRESHOLD_=65536 MALLOC_TRIM_THRESHOLD_=65536` | `sierra_nevada_*.png` (median, MAD, anomalies, spring temperature, scatter, legend) and the assembled `sierra_nevada_worked_example.png` |
| `climate/` | `pixelwise_climate_correlations` | the ERA5-Land store and the public pyramid (Equal Earth derived on the fly), GMBA (web) | a correlations zarr under `scratch/` (local, exploratory) |
| `climate/temperature_sensitivity_comparison/` | `evan_eisenman_2021_comparison` | the SNOTEL/CCSS daily archive and a WeatherBench2 ERA5 climatology cache (both self-downloaded into its `data/`), the mountain-range metrics table | `ee2021_*.png` (the recreated figures of Evan & Eisenman 2021), `runoff_onset_sensitivity_vs_ee2021.png`, `evan_eisenman_2021_snotel_station_metrics.csv` |

Order that matters: each unit's `0_aggregate` notebook before anything else in its folder; `river_basins/snow_water`
before `basin_onset` if the precipitation-as-snow panels should be filled; the two composite notebooks build their own
sweeps (skip-if-complete, `REBUILD_SWEEP`). `pixi run aggregate` runs the three aggregation notebooks headlessly.

## Conventions

- **Every read is a visible open.** `*_ds = xr.open_dataset(path)`, `*_df = pd.read_csv(path)`, `*_gdf = gpd.read_file(...)`
  with the path built from `gsro_analysis.paths`; the variable is displayed in the same cell so a run shows what was read.
  Names carry the object type: `_ds` (Dataset), `_da` (DataArray), `_df` (DataFrame), `_gdf` (GeoDataFrame).
- **Every computation the reader needs is in a cell**, one operation per cell with its constants named at the top. The
  package holds only the pipeline engine, the partial-sum math (`aggregate.reduce_partials`), three count-weighted
  helpers (`aggregate.weighted_mean` / `collapse` / `threshold`), the polar-panel styling (`plotting`), the annotated
  colorbars and the two world-map renderers.
- The mountain-range notebooks share one "analyses' view" block: CHILI classes folded, bins with ≤ 100 pixels masked, a
  bin-year needs > 30 % of the bin's median pixels, the tropical-Andes rule, and a range-year mean anomaly needs ≥ 10 %
  of the range's median pixels (the same rule as the metrics table, so the insets, the choropleths and the tables agree).
- Every quoted number exists in a `results/<version>/` table stamped with `_version`, `_git_sha` (production package),
  `_analysis_git_sha` (this repo) and `_written_at` (`gsro_analysis.results.provenance`).
- Figures are PNG at ≤ 300 dpi, saved with an explicit `savefig` path; no `plt.show()`. Regressions on the eleven
  water years report the Theil–Sen slope next to OLS (see the analysis note in `docs/aggregation_lineage.md`).
- The version comes from `settings.load_config()` only; every artifact names its version and its pixel filter
  (`fcf_lte_50` everywhere, the case study included).

Figures: one tracked set per dataset version under `analyses/<unit>/figures/<version>/` (`v9` is the frozen set from the
2025 workflow, minus the renders no code produces any more; `v10` follows the campaign); the sweeps are gitignored and
regenerable. `GSRO_OUTPUT_ROOT` redirects everything a notebook writes (figures, results, `data/aggregation`) and
`GSRO_PARTIALS_ROOT` the partials cache, so a notebook can run against a test set without touching the tree.

Credentials: none for the topic notebooks of `continents/`, `mountain_ranges/` and `river_basins/` (cubes and tables are
local) and for the Evan & Eisenman comparison; the Azure SAS token for the aggregation notebooks (partials download and
ERA5-Land zonal means; they run without it on the local cache) and for `case_studies/` and `climate/`; Earth Engine for
`river_basins/snow_water` (`settings.initialize_earthengine()`).
