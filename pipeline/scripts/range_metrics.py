"""Stage 3 (local, seconds): the per-range metrics as ONE provenance-stamped
results table.

Reads aggregated_results/<version>/mountain_ranges/all_mountain_ranges_<filter>.nc
(with the ERA5 zonal variables merged by reduce_partials.py) and writes

  analyses/mountain_ranges/results/<version>/mountain_range_metrics.csv
      one row per range: GMBA id, continent, centroid, pixel count, mean MAD and
      median onset, lapse rates (both definitions in gsro_analysis.stats),
      the mean onset anomaly per water year, and — when the ERA5 zonal means
      are in the cube — the spring-temperature sensitivity (OLS + Theil-Sen).
      Provenance columns: _git_sha (production package), _analysis_git_sha
      (this repo), _written_at.

The notebooks under analyses/ only PLOT this table; the world maps and the
per-range comparison join it to the GMBA polygons on read
(stats.range_metrics_gdf). The three legacy-named subsets and the stats
geojson that used to be written alongside were dropped on 2026-09-03.
"""

import argparse

import pandas as pd

from gsro_analysis import aggregate, paths, settings, stats
from gsro_analysis.results import save_result_table   # stamps _analysis_git_sha too


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--filter', default='fcf_lte_50')
    p.add_argument('--min-pixels', type=int, default=100,
                   help='bins with fewer pixels are masked before any statistic (default 100)')
    p.add_argument('--min-year-fraction', type=float, default=0.1,
                   help="a range-year needs this fraction of the range's median-pixel count (default 0.1)")
    return p.parse_args()


def main():
    args = parse_args()
    config = settings.load_config(args.config)
    ds = aggregate.open_aggregate('mountain_ranges', config.version, args.filter)
    ds = aggregate.threshold(ds, args.min_pixels)
    water_years = [int(y) for y in ds.water_year.values]

    dims = ['elevation', 'aspect', 'chili_class']
    metrics = pd.DataFrame({
        'name': ds.mountain_range.values,
        'GMBA_V2_ID': ds['GMBA_V2_ID'].values,
        'continent': ds['continent'].values,
        'centroid_latitude': ds['centroid_latitude'].values,
        'centroid_longitude': ds['centroid_longitude'].values,
        'total_pixels_in_range': ds['runoff_onset_median_n'].sum(dims).values,
        'mean_mad_days': aggregate.weighted_mean(ds, 'runoff_onset_mad', dims).values,
        'mean_median_onset_dowy': aggregate.weighted_mean(ds, 'runoff_onset_median', dims).values,
    }).set_index('name')

    lapse_bins = stats.lapse_rate_weighted_bins(ds).set_index('name')
    metrics['lapse_rate_weighted_bins_per_100m'] = lapse_bins['lapse_rate']
    metrics['lapse_rate_weighted_bins_r2'] = lapse_bins['r_squared']
    prof = stats.lapse_rate_profile(ds)
    metrics['snowmelt_lapse_rate_per_100m'] = prof['lapse_rate_per_100m']
    metrics['snowmelt_lapse_rate_corr'] = prof['correlation']
    metrics['snowmelt_lapse_rate_n'] = prof['n']

    anom = stats.range_mean_anomaly(ds, min_year_fraction=args.min_year_fraction)
    for y in water_years:
        metrics[f'runoff_onset_anomaly_WY{y}'] = anom.sel(water_year=y).values

    if 'temperature_2m' in ds:
        metrics = metrics.join(stats.spring_temperature_sensitivity(ds, min_year_fraction=args.min_year_fraction))
    else:
        print("note: no ERA5 variables in the cube (run era5_zonal.py, then reduce_partials.py) "
              "- temperature sensitivity skipped")

    results_dir = paths.resultsdir('mountain_ranges', config.version)
    out = save_result_table(metrics.reset_index().round(3), 'mountain_range_metrics', results_dir=results_dir)
    print(f"{len(metrics)} ranges, {len(metrics.columns)} metrics -> {out}")


if __name__ == '__main__':
    main()
