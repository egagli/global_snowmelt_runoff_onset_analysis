"""'Get ancillary data' worker: for a batch of tiles, build the eight ancillary layers on the
tile's window of the DATASET grid and commit them to the grid generation's icechunk repository
(datacube.build_ancillary_window -> write_ancillary_tile). One matrix job runs one batch;
equally runnable locally (--tile / --tile-list, --local-store for a local repo).

Tiles that already have a commit are skipped (the ledger is folded once per job). Failure =
exception = no commit: the tile is re-listed by the next dispatch. The job exits nonzero if any
tile failed, so the run shows red, but one bad tile never blocks its batchmates.

Needs Earth Engine (CHILI) + an Azure SAS token; the repository must exist (the plan job's
`get_remaining_work.py --stage ancillary` creates it).
"""

import argparse
import json
import sys
import time
import traceback

from gsro_analysis import datacube, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--manifest', help='manifest json from get_remaining_work.py --stage ancillary')
    p.add_argument('--batch', type=int, help='batch index into the manifest')
    p.add_argument('--tile', nargs=2, type=int, action='append', metavar=('ROW', 'COL'))
    p.add_argument('--tile-list', help='file with one "row col" pair per line')
    p.add_argument('--ee-key', help='Earth Engine service-account key json '
                   '(default: the production clone\'s config/ee_key.json)')
    p.add_argument('--local-store', help='path of a local icechunk repository (tests)')
    return p.parse_args()


def load_tiles(args):
    tiles = [tuple(t) for t in (args.tile or [])]
    if args.tile_list:
        with open(args.tile_list) as f:
            tiles += [tuple(map(int, line.split())) for line in f if line.strip() and not line.startswith('#')]
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
    config = settings.load_config(args.config)
    log = lambda m: print(f"    {m}", flush=True)  # noqa: E731

    print(f"Earth Engine initialized as {settings.initialize_earthengine(key_file=args.ee_key)}", flush=True)
    # warm the vector-source cache ONCE before looping: per-tile reads must never hit the network
    settings.cached_source(settings.GMBA_URL)
    settings.cached_source(settings.CONTINENTS_URL)
    settings.cached_source(settings.BASIN_ATLAS_URL, filename='BasinATLAS_Data_v10.gdb.zip',
                           expected_md5=settings.BASIN_ATLAS_MD5)

    repo = datacube.open_ancillary_repo(config, args.local_store)
    done = datacube.completed_ancillary_tiles(config, repo)
    failed = []
    for i, (row, col) in enumerate(tiles):
        t0 = time.time()
        try:
            if (row, col) in done:
                print(f"[{i+1}/{len(tiles)}] tile {row},{col}: ancillary committed, skip", flush=True)
                continue
            print(f"[{i+1}/{len(tiles)}] tile {row},{col}: start", flush=True)
            ds = datacube.build_ancillary_window(config, row, col, log=log)
            datacube.write_ancillary_tile(config, repo, row, col, ds, duration_s=time.time() - t0, log=log)
            print(f"[{i+1}/{len(tiles)}] tile {row},{col}: done ({time.time() - t0:.0f}s)", flush=True)
        except Exception:
            print(f"[{i+1}/{len(tiles)}] tile {row},{col}: FAILED after {time.time() - t0:.0f}s", flush=True)
            traceback.print_exc()
            failed.append((row, col))

    print(f"batch done: {len(tiles) - len(failed)} ok, {len(failed)} failed", flush=True)
    if failed:
        sys.exit(f"failed tiles: {failed}")


if __name__ == '__main__':
    main()
