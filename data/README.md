# data — where every input lives and how it is rebuilt

Inputs live next to the analysis that uses them; this folder holds only the one basemap every map shares. Everything
below is **input**: either downloaded by the code on first use or built by a script in `pipeline/scripts/`, so a fresh
clone can recreate all of it. Pipeline and analysis **outputs** live in `partials/<version>/` (the fleet product, cached)
and `analyses/<unit>/data/aggregation/<version>/`, `results/<version>/`, `figures/<version>/`. Sizes as of September 2026.

| Path | Size | Tracked | Provenance and how to rebuild |
| --- | --- | --- | --- |
| `data/global_hillshade_robinson.tif` | 252 MB | no | Natural Earth *Gray Earth with Shaded Relief, Hypsography, Ocean Bottom, and Drainages* (10 m raster, public domain, `settings.HILLSHADE_URL`) reprojected to World Robinson (ESRI:54030) at 1 km with gdalwarp. `pixi run hillshade` (`pipeline/scripts/get_hillshade.py`) downloads the zip into `data/sources/` and builds it; no credentials; byte-identical to the file used before the script existed. Basemap of the two inset world maps (`gsro_analysis.world_maps`), the basin maps and the case study. |
| `data/sources/` | 350 MB | no | The Natural Earth zip and its extracted GeoTIFF (the hillshade's source). |
| `analyses/continents/data/gtopo30_lat_elev_histogram.nc` | 40 KB | **yes** | Land-pixel count per continent × 1° latitude × 100 m elevation from GTOPO30 on Earth Engine, merged into the continents cube as `dem_pixel_count`. `pixi run gtopo30` (`pipeline/scripts/get_gtopo30_histogram.py`, Earth Engine key) rebuilds it; the 2026-09-04 rebuild reproduced the tracked file cell for cell. Grid- and version-independent. |
| `analyses/continents/data/geometries/continents.zip` | 1 MB | no | USGS continents (`settings.CONTINENTS_URL`), the fleet's cache (`settings.continents_zip()`); the notebooks read the same URL straight from the web. |
| `analyses/mountain_ranges/data/geometries/GMBA_Inventory_v2.0_standard_300.zip` | 33 MB | no | GMBA Mountain Inventory v2.0, standard 300 (`settings.GMBA_URL`), the fleet's cache (`settings.gmba_zip()`); the notebooks read the same URL straight from the web. |
| `analyses/mountain_ranges/label_layout.csv` | 5 KB | **yes** | The curated label set and hand-placed anchors of the two inset world maps (edited, never generated). |
| `analyses/river_basins/data/geometries/BasinATLAS_Data_v10.gdb.zip` | 2.7 GB | no | BasinATLAS v1.0 (`settings.BASIN_ATLAS_URL`, md5-checked), the one vector source that must be cached: figshare sits behind a bot challenge and GDAL cannot stream the gdb. `settings.basin_atlas_gdb()` downloads it on first use; the fleet stores HydroBASINS level 6, the analyses read level 5. The fleet jobs and the CI restore it from the Actions cache. |
| `analyses/river_basins/data/geometries/hydrobasins_level6_population.csv` | 0.5 MB | **yes** | GPW v4.11 population count (CIESIN, 2020, 30 arc-second) summed per HydroBASINS v1 level-6 basin on Earth Engine. `pixi run population` (`pipeline/scripts/get_basin_population.py`, Earth Engine key) rebuilds it; level 5 is the exact sum by `PFAF_ID // 10` (the river-basin aggregation notebook does it). Replaces the 298 MB level-5 geojson exported from the Code Editor in April 2025. |
| `analyses/climate/temperature_sensitivity_comparison/data/snotel_ccss_archive/` | 462 MB | no | Daily SNOTEL/CCSS station archive from [egagli/snotel_ccss_stations](https://github.com/egagli/snotel_ccss_stations), self-downloaded by the Evan & Eisenman notebook on first run. |
| `analyses/climate/temperature_sensitivity_comparison/data/evan_eisenman_2021_era5_T0_T1_climatology.nc` | 1 MB | no | The WeatherBench2 ERA5 (1990–2019) annual-cycle climatology fit, built by the same notebook on first run from a ~770 MB anonymous read. |
| `pipeline/tile_data/ancillary_tiles_v10.txt` | 40 KB | **yes** | The fleet work list: every tile whose composites hold data, from the dataset store's commit history (regenerated in `pipeline/pipeline.ipynb`). |

Read live, never cached: the World Bank *Major Rivers of the World* (`settings.MAJOR_RIVERS_URL`, the rivers overlay of
the basin maps) and the NOAA PSL Niño 3.4 / PDO / PNA monthly indices (`mountain_ranges/teleconnections.ipynb`).

On Azure (`snowmelt/snowmelt_runoff_onset_analysis/`, SAS token): the ERA5-Land icechunk repository per version, the
ancillary icechunk repository per grid and the per-tile partials, all built by the three GitHub Actions workflows
(README). The runoff-onset dataset store itself belongs to the production repo.
