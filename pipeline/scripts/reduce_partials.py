"""Stage 2 (local, minutes): reduce the fleet's per-tile partial sums into the
aggregate cubes the analyses read.

  aggregated_results/<version>/<group>/all_<group>_<filter>.nc
    group: by default mountain_ranges, river_basins (HydroBASINS level 5,
    derived from the stored level-6 ids) and continents (latitude x
    elevation); on demand --groups river_basins_l6 and/or continents_aspect
    (aggregate.GROUPS). One file per filter tag (aggregate.FILTERS: fcf_lte_50
    and full_dataset), compressed netCDF.

Inputs: snowmelt/snowmelt_runoff_onset_analysis/partials/<version>/tile_*.parquet (downloaded once
into aggregated_results/<version>/partials/ — re-runs only fetch new tiles),
GMBA + USGS continents (cached vector sources), the static GTOPO30 histogram
(data/gtopo30_lat_elev_histogram.nc, --build-gtopo30 regenerates it via
Earth Engine) and, if present, the ERA5 zonal files written by era5_zonal.py
(merged into the mountain-range cube). --mirror uploads the finished cubes to
snowmelt/snowmelt_runoff_onset_analysis/aggregated/<version>/ so they exist somewhere other than
this machine. Needs an Azure SAS token.
"""

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import xarray as xr

from gsro_analysis import aggregate, paths, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--filters', default='',
                   help='comma-separated filter tags (default: every aggregate.FILTERS tag present '
                        'in the partials)')
    p.add_argument('--groups', default=','.join(aggregate.DEFAULT_GROUPS),
                   help=f'comma-separated output groups from {list(aggregate.GROUPS)} '
                        f'(default: {",".join(aggregate.DEFAULT_GROUPS)})')
    p.add_argument('--build-gtopo30', action='store_true',
                   help='(re)build data/gtopo30_lat_elev_histogram.nc via Earth Engine')
    p.add_argument('--mirror', action='store_true',
                   help='upload the cubes to snowmelt/snowmelt_runoff_onset_analysis/aggregated/<version>/')
    p.add_argument('--no-download', action='store_true',
                   help='use the local partials cache as is (no Azure listing)')
    return p.parse_args()


def sync_partials(config, cache_dir, workers=16):
    """Download partials parquets missing from the local cache; drop cached
    files the ledger no longer lists (a redone tile). Returns local paths."""
    fs = config.azure_blob_fs
    prefix = f"{settings.PARTIALS_PREFIX}/{config.version}"
    remote = {p.rsplit('/', 1)[-1]: p for p in fs.ls(prefix, detail=False)
              if p.endswith('.parquet')}
    local = {f for f in os.listdir(cache_dir) if f.endswith('.parquet')}
    stale = local - set(remote)
    for f in stale:
        os.remove(cache_dir / f)
    missing = sorted(set(remote) - local)
    print(f"partials on Azure: {len(remote)} | cached: {len(local) - len(stale)} | "
          f"fetching {len(missing)} | dropped stale {len(stale)}", flush=True)

    def fetch(name):
        fs.get_file(remote[name], str(cache_dir / f"{name}.part"))
        os.replace(cache_dir / f"{name}.part", cache_dir / name)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(fetch, missing))
    return sorted(cache_dir / f for f in remote)


# the partials' key columns: every row is one tile's contribution to one cube cell
KEY_COLS = ['filter_tag', 'unit_type', 'unit_id', 'elevation', 'aspect', 'latitude', 'chili_class']


def load_partials_summed(files, batch=200):
    """Read the tiles' partials and SUM them over tiles as they are read (groupby-sum per
    batch of files, then over the batches). The reduce is a sum over identical keys, so the
    result is the same as concatenating every row first, at a fraction of the memory: the
    4,320 v10 tiles hold ~14 M rows x 85 float64 columns (~10 GB concatenated; the 15 GB
    dev box was OOM-killed on 2026-09-04) but only as many summed rows as there are populated
    cube cells. ``min_count=1`` keeps a column NaN when no tile reported it (the same rule
    aggregate.reduce_partials applies). Returns (summed frame, number of raw rows)."""
    acc, n_rows = None, 0
    for i in range(0, len(files), batch):
        chunk = pd.concat((pd.read_parquet(f) for f in files[i:i + batch]), ignore_index=True)
        n_rows += len(chunk)
        chunk = chunk.drop(columns=[c for c in ('tile_row', 'tile_col') if c in chunk])
        summed = chunk.groupby(KEY_COLS, sort=False, dropna=False).sum(min_count=1)
        acc = summed if acc is None else pd.concat([acc, summed]).groupby(level=KEY_COLS, sort=False, dropna=False).sum(min_count=1)
        del chunk, summed
    return acc.reset_index(), n_rows


def encoding_for(ds):
    enc = {}
    for v in ds.data_vars:
        e = {'zlib': True, 'complevel': 4}
        if ds[v].dtype.kind == 'f':
            e['dtype'] = 'float32'
        enc[v] = e
    return enc


def main():
    args = parse_args()
    config = settings.load_config(args.config)
    version = config.version
    t0 = time.time()

    cache_dir = paths.partials_cache(version)
    if args.no_download:
        files = sorted(cache_dir / f for f in os.listdir(cache_dir) if f.endswith('.parquet'))
    else:
        files = sync_partials(config, cache_dir)
    if not files:
        sys.exit(f"no partials for {version} under {settings.PARTIALS_PREFIX}/{version}")
    partials, n_rows = load_partials_summed(files)
    n_tiles = len(files)
    print(f"summed {n_rows:,} partial rows from {n_tiles} tiles into {len(partials):,} cube cells "
          f"({time.time() - t0:.0f}s)", flush=True)

    present = set(partials['filter_tag'].unique())
    filters = [f for f in args.filters.split(',') if f] or [t for t in aggregate.FILTERS if t in present]
    groups = [g for g in args.groups.split(',') if g]
    unknown = set(groups) - set(aggregate.GROUPS)
    if unknown:
        sys.exit(f"unknown groups {sorted(unknown)}; choose from {list(aggregate.GROUPS)}")
    water_years = [int(y) for y in config.water_years]

    gmba_gdf = aggregate.load_gmba()
    continents_gdf = aggregate.load_continents()
    gtopo_path = paths.gtopo30_histogram()
    if args.build_gtopo30 or not gtopo_path.exists():
        print("building the GTOPO30 latitude x elevation histogram (Earth Engine) ...", flush=True)
        aggregate.gtopo30_lat_elev_histogram(continents_gdf).to_netcdf(
            gtopo_path, encoding={'pixel_count': {'zlib': True, 'complevel': 4}})
    gtopo30 = xr.open_dataset(gtopo_path)

    zonal = {}
    for group in groups:
        if group.startswith('continents'):
            continue
        zp = paths.era5_zonal(version, group)
        if zp.exists():
            zonal[group] = xr.open_dataset(zp)
        else:
            print(f"note: no ERA5 zonal file for {group} ({zp}); run era5_zonal.py --units {group} "
                  f"to add the climate variables", flush=True)

    written = []
    for group in groups:
        for tag in filters:
            t1 = time.time()
            try:
                ds = aggregate.reduce_partials(partials, group, tag, water_years)
            except ValueError as e:
                print(f"  {group}/{tag}: {e}", flush=True)
                continue
            if group == 'mountain_ranges':
                ds = aggregate.finalize_mountain_ranges(ds, gmba_gdf, continents_gdf,
                                                       era5_zonal_ds=zonal.get(group))
            elif group.startswith('continents'):
                ds = aggregate.finalize_continents(ds, gtopo30)
            elif group.startswith('river_basins') and group in zonal:
                z = zonal[group].rename({'PFAF_ID': 'river_basin'})
                z = z.sel(river_basin=z['river_basin'].isin(ds['river_basin'].values))
                ds = xr.merge([ds, z], combine_attrs='drop_conflicts', join='left')
            ds.attrs.update({'dataset_version': version, 'n_tiles': int(n_tiles),
                             'produced_by': 'pipeline/scripts/reduce_partials.py'})
            out = paths.aggregate(group, version, tag)
            out.parent.mkdir(parents=True, exist_ok=True)
            tmp = out.with_suffix('.nc.tmp')
            ds.to_netcdf(tmp, encoding=encoding_for(ds))
            os.replace(tmp, out)
            written.append(out)
            print(f"  wrote {out.relative_to(paths.ROOT) if out.is_relative_to(paths.ROOT) else out} "
                  f"{os.path.getsize(out) / 1e6:.1f} MB, dims {dict(ds.sizes)} ({time.time() - t1:.0f}s)",
                  flush=True)

    if args.mirror and written:
        fs = config.azure_blob_fs
        for out in written:
            dest = f"{settings.AGGREGATED_PREFIX}/{version}/{out.parent.name}/{out.name}"
            fs.put_file(str(out), dest)
            print(f"  mirrored -> {dest}", flush=True)
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
