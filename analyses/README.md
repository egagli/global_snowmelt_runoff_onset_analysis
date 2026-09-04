# analyses

One folder per aggregation unit (each reads the cube of that unit type, written by
`pipeline/scripts/reduce_partials.py`), plus the Sierra Nevada case study (reads the dataset store
directly) and `climate/` (ERA5-Land products and a comparison with a published sensitivity model).
Every notebook opens with a markdown cell that says what it reads and writes; run its setup cells
(imports, `config = settings.load_config()`, the cube) first, then any section. Display-range lists
are filtered to what exists in the cube, so a partially processed dataset version runs.

| Folder | Notebook | Reads | Writes (`figures/<version>/`, `results/<version>/`) |
| --- | --- | --- | --- |
| `global/` | `lat_elev_binning` | continents cube (`fcf_lte_50` and, for the share of pixels above 5000 m, `full_dataset`), `data/gtopo30_lat_elev_histogram.nc` | `continent_lat_and_elev_binned.png` (median, MAD, sunny–shaded difference per continent, annotated colorbars drawn in place), `continent_lat_and_elev_binned_anomaly.png`, FCF/CHILI correlation panels |
| `global/` | `sunny_shaded` | continents cube | `global_scatter_latitude_vs_onset.png` (seasonal modulation of the sunny–shaded offset), offset-vs-onset and MAD-vs-onset scatters |
| `global/` | `lapse_rates` | continents cube | `lapse_rates_by_latitude_groupby_continents.png`, `global_lapse_rates_by_continent.png`; latitude + elevation regressions |
| `mountain_ranges/` | `topography` | mountain-range cube | `polar_triplets_all_ranges.png`, `global_lapse_rates_by_mountain_range.png`, `global_mad_by_mountain_range.png` |
| `mountain_ranges/` | `topography_triplet_composite_figure` | mountain-range cube, `results/<version>/mountain_range_metrics.csv` (joined to the GMBA polygons on read), `label_layout.csv`, `data/global_hillshade_robinson.tif` | the per-range triplet sweep `triplets/pngs/` (gitignored), `polar_triplet_legend.png`, the composite `global_lapse_rates_with_triplets_map.png` (`gsro_analysis.world_maps`) |
| `mountain_ranges/` | `anomalies` | mountain-range cube, GMBA polygons | `snowmelt_onset_anomalies_by_mountain_range.png`, `regional_panels/*.png` (gitignored sweep), `regional_panels/combined_annual/*.png`, `median_runoff_and_anomaly_subset.png` |
| `mountain_ranges/` | `temperature_sensitivity` | mountain-range cube (with the ERA5-Land zonal means), `results/<version>/mountain_range_metrics.csv` | `snowmelt_onset_spring_temp_sensitivity_histograms.png`, `mountain_range_spring_temp_vs_runoff_onset_anomaly_scatterplots_theil_sen.png`, regression heatmaps |
| `mountain_ranges/` | `spring_temperature_sensitivity_composite_figure` | as the triplet composite, plus the ERA5-Land zonal means | the per-range scatterplot sweep `anomaly_scatterplots/pngs/` (gitignored), the composite `global_temperature_sensitivity_map.png` |
| `mountain_ranges/` | `teleconnections` | NOAA PSL ENSO / PDO / PNA indices (online) | `teleconnection_indices.png` |
| `river_basins/` | `basin_onset` | river-basin cube (HydroBASINS level 5), BasinATLAS level-5 polygons (cached gdb), `data/geometries/Hydrobasins_L5_Population_Global.geojson`, the hillshade, `results/<version>/river_basin_snow_water.csv` if present, the rivers overlay if present | `global_median_runoff_onset_and_mad.png`, `global_runoff_onset_anomaly_by_basin.png`, `{hma,wus}_*` regional panels |
| `river_basins/` | `snow_water` | river-basin cube, ERA5-Land via Earth Engine (April-1 / October-1 SWE, precipitation-as-snow share) | `results/<version>/river_basin_snow_water.csv` (read by `basin_onset`), SWE / water-volume vs onset scatters |
| `case_studies/sierra_nevada/` | `sierra_nevada` | the dataset store (4x coarsened), the ERA5-Land store, the hillshade, forest cover (FCF ≤ 50 % masked like the pipeline) | `sierra_nevada_*.png` (median, MAD, anomalies, spring temperature, scatter, legend) and the assembled `sierra_nevada_worked_example.png` |
| `climate/` | `pixelwise_climate_correlations` | ERA5-Land store and the public pyramid (Equal Earth derived on the fly) | a correlations zarr under `scratch/` (local, exploratory) |
| `climate/temperature_sensitivity_comparison/` | `evan_eisenman_2021_comparison` | SNOTEL/CCSS daily archive (self-downloaded into `data/snotel_ccss_archive/`), an ERA5 T₀/T₁ climatology cache (`data/era5/`), `results/<version>/mountain_range_metrics.csv` | `ee2021_*.png` (the recreated figures of Evan & Eisenman 2021), `runoff_onset_sensitivity_vs_ee2021.png`, `results/<version>/evan_eisenman_2021_snotel_station_metrics.csv` |

Order that matters: `pipeline/scripts/range_metrics.py` before anything that reads the metrics table;
`river_basins/snow_water` before `basin_onset` if the precipitation-as-snow panels should be filled;
the two composite notebooks build their own sweeps (skip-if-complete, `REBUILD_SWEEP`).

Conventions: the version comes from `settings.load_config()` only; every artifact names its version
and its pixel filter (`fcf_lte_50` everywhere, the case study included); statistics come from the cube
via `aggregate.weighted_mean` / `collapse` / `threshold`, never re-derived by hand; regressions on the
eleven water years report the Theil–Sen slope next to OLS (see the analysis note in `docs/aggregation_lineage.md`); the mountain-range
notebooks start from `stats.prepare_mountain_ranges` (CHILI collapsed, bins with ≤ 100 pixels masked,
a bin-year needs > 30 % of the bin's median pixels, the tropical-Andes rule); renders are PNG at
≤ 300 dpi; annotated colorbars come from `gsro_analysis.colorbars`; the two inset world maps come from
`gsro_analysis.world_maps` and `mountain_ranges/label_layout.csv` (the curated label set and anchors,
the one figure input that is neither data nor code).

Figures: one tracked set per dataset version under `analyses/<unit>/figures/<version>/` (`v9` is the
frozen set from the 2025 workflow; `v10` follows the campaign); the sweeps are gitignored and
regenerable. Results tables carry `_git_sha` (production package), `_analysis_git_sha` (this repo) and
`_written_at`. `GSRO_OUTPUT_ROOT` redirects figures and results, `GSRO_AGGREGATED_ROOT` the cubes,
so a notebook can run against a test set without touching the tree.

Credentials: none for `global/`, `mountain_ranges/` (cubes and tables are local) and the Evan &
Eisenman comparison; the Azure SAS token for the store and ERA5 stack (`case_studies/`, `climate/`);
Earth Engine for `river_basins/snow_water` (`settings.initialize_earthengine()`).
