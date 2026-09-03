"""Stage 2b: per-unit ERA5-Land anomaly zonal means (mountain ranges AND
river basins) from the version's anomaly store — the climate side of the
mountain-range cube (spring-temperature sensitivity) and of the river-basin
cube. Unit types: mountain_ranges (GMBA), river_basins (HydroBASINS level 5,
the default basin cube) and river_basins_l6 (level 6, for the on-demand
level-6 cube); the polygons come from the cached BasinATLAS gdb.

Weight per ERA5 cell = polygon coverage fraction x seasonal-snow fraction x
onset validity (era5.zonal_anomalies); one streamed pass over the store per
unit type (~1 GB peak, minutes — runs on a laptop; no Earth Engine). Writes
  aggregated_results/<version>/era5_zonal/era5_anomaly_<unit_type>.nc
which reduce_partials.py merges into the cubes. Needs an Azure SAS token.
"""

import argparse
import time

from gsro_analysis import aggregate, era5, paths, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--units', default='mountain_ranges,river_basins',
                   help='comma-separated: mountain_ranges, river_basins (level 5), river_basins_l6')
    p.add_argument('--variables', default='',
                   help='comma-separated ERA5 variables (default: all 8)')
    p.add_argument('--overwrite', action='store_true')
    return p.parse_args()


def load_units(unit_type):
    if unit_type == 'mountain_ranges':
        return aggregate.load_gmba(), 'GMBA_V2_ID', aggregate.RANGE_GEOMETRY_FIXES, aggregate.SKIP_RANGES
    if unit_type in ('river_basins', 'river_basins_l6'):
        import geopandas as gpd
        level = 5 if unit_type == 'river_basins' else 6
        basins = gpd.read_file(settings.basin_atlas_gdb(), layer=settings.basin_atlas_layer(level))
        # PFAF_ID is unique per basin except for a handful of multipart
        # records: dissolve so one row = one id (as the ancillary PFAF_ID
        # layer treats them)
        if basins['PFAF_ID'].duplicated().any():
            basins = basins.dissolve(by='PFAF_ID', as_index=False)
        return basins, 'PFAF_ID', None, ()
    raise SystemExit(f"unknown unit type {unit_type}")


def main():
    args = parse_args()
    config = settings.load_config(args.config)
    variables = [v for v in args.variables.split(',') if v] or None
    anomaly = era5.open_anomaly(config)
    for unit_type in [u for u in args.units.split(',') if u]:
        out = paths.era5_zonal(config.version, unit_type)
        if out.exists() and not args.overwrite:
            print(f"skip {unit_type} (exists: {out})")
            continue
        t0 = time.time()
        units_gdf, id_col, fixes, skip = load_units(unit_type)
        print(f"{unit_type}: {len(units_gdf)} polygons", flush=True)
        ds = era5.zonal_anomalies(config, units_gdf, id_col, variables=variables,
                                  geometry_fixes=fixes, skip_names=skip,
                                  anomaly_ds=anomaly,
                                  progress=lambda m: print(m, flush=True))
        ds.attrs['unit_type'] = unit_type
        tmp = out.with_suffix('.nc.tmp')
        ds.to_netcdf(tmp, encoding={v: {'zlib': True, 'complevel': 4} for v in ds.data_vars})
        tmp.replace(out)
        print(f"wrote {out} dims {dict(ds.sizes)} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == '__main__':
    main()
