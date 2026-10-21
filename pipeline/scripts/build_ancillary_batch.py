"""'Get ancillary data' worker: for a batch of tiles, build the ancillary layers on the tile's
window of the DATASET grid and commit them to the grid generation's icechunk repository. One
matrix job runs one batch; equally runnable locally (--tile / --tile-list, --local-store for a
local repo).

Per tile, the ledger fold (once per job) says which of datacube.ANCILLARY_LAYERS are committed:
  all of them                      -> skip
  none                             -> full build (build_ancillary_window -> write_ancillary_tile)
  all but some unit-id layers      -> refresh only those (refresh_unit_layers: polygons, no Earth
                                      Engine) — how a layer added after the campaign (GRDC_MRBID,
                                      2026-10-09) is backfilled
  anything else                    -> full build
Failure = exception = no commit: the tile is re-listed by the next dispatch. The job exits nonzero
if any tile failed, so the run shows red, but one bad tile never blocks its batchmates.

Needs an Azure SAS token, and Earth Engine (CHILI) for full builds only; the repository must exist
(the plan job's `get_remaining_work.py --stage ancillary` creates it and adds missing layers).
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

    # warm the vector-source cache ONCE before looping: per-tile reads must never hit the network
    settings.gmba_zip()
    settings.continents_zip()
    settings.basin_atlas_gdb()
    settings.grdc_mrb_zip()
    settings.water_towers_zip()          # 990 MB, schema 3 (WTU_SIGNED)

    repo = datacube.open_ancillary_repo(config, args.local_store)
    committed = datacube.committed_ancillary_layers(repo)
    ee_ready = False
    failed = []
    for i, (row, col) in enumerate(tiles):
        t0 = time.time()
        try:
            have = committed.get((row, col), set())
            missing = tuple(n for n in datacube.ANCILLARY_LAYERS if n not in have)
            if not missing:
                print(f"[{i+1}/{len(tiles)}] tile {row},{col}: ancillary committed, skip", flush=True)
                continue
            if have and set(missing) <= set(datacube.UNIT_ID_LAYERS):
                print(f"[{i+1}/{len(tiles)}] tile {row},{col}: refresh {missing} (no Earth Engine)", flush=True)
                datacube.refresh_unit_layers(config, row, col, layers=missing, repo=repo, log=log)
            else:
                if not ee_ready:
                    print(f"Earth Engine initialized as {settings.initialize_earthengine(key_file=args.ee_key)}", flush=True)
                    ee_ready = True
                print(f"[{i+1}/{len(tiles)}] tile {row},{col}: start (full build)", flush=True)
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
