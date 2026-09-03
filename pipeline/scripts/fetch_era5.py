"""Fetch ERA5-Land water-year stores (stage: analyses/climate step 0).

One store per water year under era5_data/<version>/ — skip-if-exists, so
listing several years only fetches the missing ones. Per-variable writes
keep peak memory ~1/8 of a whole year; still, the fetch wants a machine
with >= 8 GB free — the 16 GB GitHub runners qualify, the WSL dev box does
not (.github/workflows/era5_acquire.yml dispatches this).

For a NEW dataset version, copy the previous version's year stores first
(they are version-independent acquisitions) and fetch only the new year:
    python fetch_era5.py --copy-from-version v10 --water-years 2026

Needs Earth Engine + an Azure SAS token.
"""

import argparse

from global_snowmelt_runoff_onset.config import Config

from gsro_analysis import era5, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--water-years', default='',
                   help='comma-separated, e.g. "2025" or "2024,2025"; '
                        'default = every year in config.water_years')
    p.add_argument('--copy-from-version', default='',
                   help='server-side copy this version\'s existing year '
                        'stores into the target prefix before fetching '
                        '(skips years that already exist)')
    p.add_argument('--ee-key', help='Earth Engine service-account key json '
                   '(default: the production clone\'s config/ee_key.json)')
    p.add_argument('--build-anomaly', action='store_true',
                   help='after the year stores, (re)build the anomaly store '
                        '(all water years as the base period; ~4 GB)')
    return p.parse_args()


def copy_stores(config, from_version):
    """Server-side copy era5_water_year_*.zarr from another version's
    prefix into this version's prefix."""
    fs = config.azure_blob_fs
    src_prefix = f"{settings.ERA5_DATA_PREFIX}/{from_version}"
    for src in sorted(fs.glob(f"{src_prefix}/era5_water_year_*.zarr")):
        year = int(src.rsplit('_', 1)[-1].removesuffix('.zarr'))
        dst = era5.water_year_store_path(config, year)
        if era5.verify_and_mark(config, dst):
            print(f"skip copy WY{year} (complete: {dst})")
            continue
        print(f"copy WY{year}: {src} -> {dst}")
        for f in fs.find(src):
            fs.cp_file(f, f.replace(src, dst))
        if not era5.verify_and_mark(config, dst):
            raise RuntimeError(f"WY{year}: copy verification failed ({dst})")


def main():
    args = parse_args()
    config = Config(args.config)
    settings.initialize_earthengine(key_file=args.ee_key)
    if args.copy_from_version:
        copy_stores(config, args.copy_from_version)
    years = ([int(y) for y in args.water_years.split(',') if y]
             if args.water_years else None)
    era5.build_water_year_stores(config, water_years=years)
    if args.build_anomaly:
        era5.build_anomaly_store(config)


if __name__ == '__main__':
    main()
