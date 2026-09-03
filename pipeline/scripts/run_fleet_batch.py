"""Fleet batch worker: stages 0 (ancillary) + 1 (tabulate -> partial sums)
for a batch of tiles.

The GitHub Actions entry point (one matrix job runs one batch), equally
runnable locally. Per tile (datacube.process_tile): build the ancillary
zarr if its completion marker is missing, tabulate in memory, write the
tile's PARTIAL SUMS parquet (the aggregation input) — and, only with
--keep-pixels, the per-pixel parquet. Tiles whose partials blob exists are
skipped, so re-dispatching a batch only runs what's missing.

Failure policy (fleet rule): a tile that fails leaves NO output — no
completion marker, no partials — and the worker moves on to the next tile.
The job exits nonzero if any tile failed, so the run shows red, but one bad
tile never blocks its batchmates. Re-dispatch re-lists failed tiles.

Tiles come from a manifest written by get_remaining_work.py (--manifest +
--batch), or directly via --tile / --tile-list for local runs.

Needs Earth Engine (CHILI) + an Azure SAS token (AZURE_STORAGE_SAS_TOKEN or
the production clone's config/sas_token.txt).
"""

import argparse
import json
import sys
import time
import traceback

from global_snowmelt_runoff_onset.config import Config

from gsro_analysis import datacube, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--manifest', help='manifest json from get_remaining_work.py')
    p.add_argument('--batch', type=int, help='batch index into the manifest')
    p.add_argument('--tile', nargs=2, type=int, action='append',
                   metavar=('ROW', 'COL'))
    p.add_argument('--tile-list', help='file with one "row col" pair per line')
    p.add_argument('--ee-key', help='Earth Engine service-account key json '
                   '(default: the production clone\'s config/ee_key.json)')
    p.add_argument('--keep-pixels', action='store_true',
                   help='also write the per-pixel parquet (opt-in product; '
                        '~35 MB/tile, ~150 GB for the grid)')
    return p.parse_args()


def load_tiles(args):
    tiles = [tuple(t) for t in (args.tile or [])]
    if args.tile_list:
        with open(args.tile_list) as f:
            tiles += [tuple(map(int, line.split())) for line in f
                      if line.strip() and not line.startswith('#')]
    if args.manifest is not None:
        if args.batch is None:
            sys.exit('--manifest requires --batch')
        with open(args.manifest) as f:
            manifest = json.load(f)
        tiles += [tuple(t) for t in manifest['batches'][args.batch]]
    if not tiles:
        sys.exit('no tiles given (--manifest/--batch, --tile, or --tile-list)')
    return tiles


def main():
    args = parse_args()
    tiles = load_tiles(args)
    config = Config(args.config)
    fs = config.azure_blob_fs

    identity = settings.initialize_earthengine(key_file=args.ee_key)
    print(f"Earth Engine initialized as {identity}", flush=True)

    # warm the vector-source cache ONCE before looping (fleet rule in
    # pipeline/README.md) — per-tile reads must never hit the network
    settings.cached_source(settings.GMBA_URL)
    settings.cached_source(settings.CONTINENTS_URL)
    settings.cached_source(settings.BASIN_ATLAS_URL,
                           filename='BasinATLAS_Data_v10.gdb.zip',
                           expected_md5=settings.BASIN_ATLAS_MD5)

    global_ds = config.open_runoff_onset_dataset(chunks=None,
                                                 mask_and_scale=True)
    failed = []
    for i, (row, col) in enumerate(tiles):
        t0 = time.time()
        try:
            if fs.exists(datacube.partials_tile_path(config, row, col)) and (
                    not args.keep_pixels
                    or fs.exists(datacube.parquet_tile_path(config, row, col))):
                print(f"[{i+1}/{len(tiles)}] tile {row},{col}: partials exist, skip",
                      flush=True)
                continue
            had_ancillary = datacube.ancillary_tile_complete(config, row, col)
            n_px, n_rows = datacube.process_tile(
                config, row, col, global_ds=global_ds, keep_pixels=args.keep_pixels)
            print(f"[{i+1}/{len(tiles)}] tile {row},{col}: "
                  f"{'ancillary reused' if had_ancillary else 'ancillary built'}, "
                  f"{n_px:,} px -> {n_rows:,} partial rows ({time.time() - t0:.0f}s)",
                  flush=True)
        except Exception:
            print(f"[{i+1}/{len(tiles)}] tile {row},{col}: FAILED after "
                  f"{time.time() - t0:.0f}s", flush=True)
            traceback.print_exc()
            failed.append((row, col))

    print(f"batch done: {len(tiles) - len(failed)} ok, {len(failed)} failed",
          flush=True)
    if failed:
        sys.exit(f"failed tiles: {failed}")


if __name__ == '__main__':
    main()
