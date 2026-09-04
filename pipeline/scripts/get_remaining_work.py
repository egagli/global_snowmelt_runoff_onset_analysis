"""Dispatcher of the two tile waves: derive the remaining work from the ledgers.

--stage ancillary   'Get ancillary data': tiles of the work list without a commit in the grid
                    generation's ancillary icechunk repository (the fold over its history). The
                    repository is created with the empty template if it does not exist.
--stage partials    'Process tiles to parquets': tiles WITH an ancillary commit and without a
                    partials parquet (blob existence, one list call). Tiles without an ancillary
                    commit are reported as blocked, never dispatched.

A tile is dispatched only once its inputs exist; the workers skip finished tiles themselves, so
a stale listing is harmless. Batches go to a manifest json shared between the plan job and the
matrix jobs of one run. Prints a summary to stderr and a compact json to stdout:
  {"batch_index": [0, 1, ...], "n_remaining": N, "n_batches": B}
With $GITHUB_OUTPUT set, also appends matrix=<json>, count=<B> there.

--start-fresh (the workflows' off-by-default checkbox, the only caller) first DELETES this
stage's products of the version on Azure: the ancillary repository (--stage ancillary) or the
partials and pixel tables (--stage partials), via datacube.reset_version.

Needs only the Azure SAS token.
"""

import argparse
import json
import os
import sys
from pathlib import Path

from gsro_analysis import datacube, settings

DEFAULT_TILE_LIST = (Path(__file__).resolve().parent.parent / 'tile_data' / 'ancillary_tiles_v10.txt')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--stage', choices=['ancillary', 'partials'], required=True)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--tile-list', default=str(DEFAULT_TILE_LIST),
                   help='full campaign work list, one "row col" per line')
    p.add_argument('--batch-size', type=int, default=36,
                   help='tiles per matrix job (256-job matrix limit: ceil(remaining/batch_size) must stay <= 256)')
    p.add_argument('--max-batches', type=int, default=0,
                   help='cap the number of batches (0 = no cap); local runs only, the workflows dispatch everything')
    p.add_argument('--start-fresh', action='store_true',
                   help="DELETE this stage's products of the version first (see the module docstring)")
    p.add_argument('--manifest-out', default='fleet_manifest.json')
    return p.parse_args()


def main():
    args = parse_args()
    with open(args.tile_list) as f:
        tiles = [tuple(map(int, line.split())) for line in f if line.strip() and not line.startswith('#')]
    config = settings.load_config(args.config)

    if args.stage == 'ancillary':
        if args.start_fresh:
            datacube.reset_version(config, what=('ancillary_grid',), confirm=True)
        repo = datacube.initialize_ancillary_store(config, log=lambda m: print(m, file=sys.stderr))
        done = datacube.completed_ancillary_tiles(config, repo)
        remaining = [t for t in tiles if t not in done]
        blocked = []
        summary = f"ancillary committed {len(done)}"
    else:
        if args.start_fresh:
            datacube.reset_version(config, what=('partials', 'pixel_tables'), confirm=True)
        done_anc = datacube.completed_ancillary_tiles(config)
        done_part = datacube.completed_partials_tiles(config)
        remaining = [t for t in tiles if t in done_anc and t not in done_part]
        blocked = [t for t in tiles if t not in done_anc]
        summary = f"ancillary committed {len(done_anc)} | partials done {len(done_part)} | blocked (no ancillary) {len(blocked)}"

    batches = [remaining[i:i + args.batch_size] for i in range(0, len(remaining), args.batch_size)]
    if args.max_batches:
        batches = batches[:args.max_batches]
    if len(batches) > 256:
        sys.exit(f"{len(batches)} batches exceeds the 256-job matrix limit; "
                 f"raise --batch-size (need >= {-(-len(remaining) // 256)})")

    manifest = {'config': args.config, 'version': config.version, 'stage': args.stage,
                'n_work_list': len(tiles), 'n_remaining': len(remaining), 'n_blocked': len(blocked),
                'batch_size': args.batch_size, 'batches': [[list(t) for t in b] for b in batches]}
    with open(args.manifest_out, 'w') as f:
        json.dump(manifest, f)

    print(f"stage {args.stage}: work list {len(tiles)} | {summary} | remaining {len(remaining)} | "
          f"{len(batches)} batches of <= {args.batch_size}", file=sys.stderr)

    result = {'batch_index': list(range(len(batches))), 'n_remaining': len(remaining), 'n_batches': len(batches)}
    print(json.dumps(result))
    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            f.write(f"matrix={json.dumps({'batch_index': result['batch_index']})}\n")
            f.write(f"count={len(batches)}\n")


if __name__ == '__main__':
    main()
