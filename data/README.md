# data/ — shared inputs

Everything in here is **input**, gitignored except this README and the 40 KB GTOPO30 histogram.
Pipeline and analysis **outputs** live in `aggregated_results/<version>/` and `analyses/<unit>/`.
Sizes as of September 2026.

| Path | Size | Provenance |
| --- | --- | --- |
| `gtopo30_lat_elev_histogram.nc` | 40 KB | **tracked.** Land-pixel count per continent × 1° latitude × 100 m elevation from GTOPO30 via Earth Engine (`aggregate.gtopo30_lat_elev_histogram`; `pipeline/scripts/reduce_partials.py --build-gtopo30` regenerates it). Merged into the continents cube as `dem_pixel_count`; grid- and version-independent. |
| `geometries/sources/` | 2.7 GB | The vector sources the ancillary build and the reduce use, downloaded once by `settings.cached_source`: GMBA Inventory v2.0 standard 300 (`settings.GMBA_URL`), the USGS continents (`settings.CONTINENTS_URL`), BasinATLAS v1.0 gdb (`settings.BASIN_ATLAS_URL`, md5-checked; the fleet stores HydroBASINS level 6 and the analyses read level 5). The fleet jobs restore this folder from the Actions cache. |
| `geometries/Hydrobasins_L5_Population_Global.geojson` | 298 MB | HydroBASINS level-5 basins with a summed population count (`HYBAS_ID`, `PFAF_ID`, `total_population`), read by `river_basins/basin_onset.ipynb` and `snow_water.ipynb` through `stats.basin_summary`. Made in the Earth Engine Code Editor in April 2025 with the script saved verbatim as `pipeline/scripts/hydrobasins_population_gee.js`: GPW v4.11 population count (CIESIN, most recent epoch, native 30 arc-second grid) summed per HydroBASINS v1 level-5 polygon (`WWF/HydroSHEDS/v1/Basins/hybas_5`) with `reduceRegions`, exported as GeoJSON. Regenerable; a Python port at level 6 (population is additive, so level 5 is an exact prefix sum) is planned. |
| `geometries/majorrivers_0_0/` | small | Major-rivers shapefile (ESRI Data & Maps). Its license does not allow redistribution and there is no download URL, so it cannot be obtained by a stranger: `river_basins/basin_onset.ipynb` skips the overlay when it is absent. The public-domain drop-in, if it is ever swapped, is Natural Earth's `ne_10m_rivers_lake_centerlines` (fetchable through `settings.cached_source`). |
| `global_hillshade_robinson.tif` | 241 MB | Global hillshade basemap (Natural Earth, 1 km, ESRI:54030), regenerable with the production repo's `visualize/data/` download notebook. Basemap of the two inset world maps (`gsro_analysis.world_maps`, grey stretch 1–231, value 0 = outside the ellipse) and of the basin and case-study maps. |
| `era5/` | 1 MB | `evan_eisenman_2021_era5_T0_T1_climatology.nc`, the WeatherBench2 (ERA5 0.25°, 1990–2019) annual-cycle climatology cache that `climate/temperature_sensitivity_comparison/evan_eisenman_2021_comparison.ipynb` builds on first run. All ERA5-Land data lives on Azure (`snowmelt/analysis/era5_data/<version>/`, via `gsro_analysis.era5`). |
| `snotel_ccss_archive/` | 462 MB | Daily SNOTEL/CCSS station archive, self-downloaded by the same notebook from [egagli/snotel_ccss_stations](https://github.com/egagli/snotel_ccss_stations). Regenerable. |

Two other inputs live elsewhere: the fleet work list `pipeline/tile_data/ancillary_tiles_v10.txt`
(every tile whose composites hold data, from the icechunk history; regenerated in
`pipeline/pipeline.ipynb`) and the curated label table of the two world maps,
`analyses/mountain_ranges/label_layout.csv` (hand-placed by design; edited, never generated).
