"""Fleet dispatcher: derive the remaining pipeline work from blob storage.

The ledger is blob existence — no database, no status file to drift:
- stage 0 done  = completion marker in the ancillary _complete/ ledger
  (flat directory, ONE list call; written only after a successful save)
- stage 1 done  = the tile's PARTIAL-SUMS parquet exists (single blob =
  atomic write); the optional per-pixel parquet is not part of the ledger

A tile is remaining if either stage is missing; the batch worker skips the
finished stage per tile. Batches are written to a manifest json that the
GitHub Actions run shares between the plan job and the matrix jobs, so every
job works from the same listing even as tiles complete mid-run.

Prints a summary to stderr and a compact json to stdout:
  {"batch_index": [0, 1, ...], "n_remaining": N, "n_batches": B}
With $GITHUB_OUTPUT set, also appends matrix=<json>, count=<B> there.

Needs only the Azure SAS token (one list call per stage prefix).
"""

import argparse
import json
import os
import sys
from pathlib import Path

from global_snowmelt_runoff_onset.config import Config

from gsro_analysis import datacube, settings

DEFAULT_TILE_LIST = (Path(__file__).resolve().parent.parent
                     / 'tile_data' / 'ancillary_tiles_v10.txt')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--tile-list', default=str(DEFAULT_TILE_LIST),
                   help='full campaign work list, one "row col" per line')
    p.add_argument('--batch-size', type=int, default=36,
                   help='tiles per matrix job (256-job matrix limit: '
                        'ceil(remaining/batch_size) must stay <= 256)')
    p.add_argument('--max-batches', type=int, default=0,
                   help='cap the number of batches this run (0 = no cap); '
                        'use e.g. 1 for a smoke test')
    p.add_argument('--manifest-out', default='fleet_manifest.json')
    return p.parse_args()


def main():
    args = parse_args()
    with open(args.tile_list) as f:
        tiles = [tuple(map(int, line.split())) for line in f
                 if line.strip() and not line.startswith('#')]

    config = Config(args.config)
    done_ancillary = datacube.completed_ancillary_tiles(config)

    done_partials = datacube.completed_partials_tiles(config)

    remaining = [t for t in tiles
                 if t not in done_ancillary or t not in done_partials]

    batches = [remaining[i:i + args.batch_size]
               for i in range(0, len(remaining), args.batch_size)]
    if args.max_batches:
        batches = batches[:args.max_batches]
    if len(batches) > 256:
        sys.exit(f"{len(batches)} batches exceeds the 256-job matrix limit; "
                 f"raise --batch-size (need >= {-(-len(remaining) // 256)})")

    manifest = {
        'config': args.config,
        'version': config.version,
        'n_work_list': len(tiles),
        'n_done_ancillary': len(done_ancillary),
        'n_done_partials': len(done_partials),
        'n_remaining': len(remaining),
        'batch_size': args.batch_size,
        'batches': [[list(t) for t in b] for b in batches],
    }
    with open(args.manifest_out, 'w') as f:
        json.dump(manifest, f)

    print(f"work list {len(tiles)} | ancillary done {len(done_ancillary)} | "
          f"partials done {len(done_partials)} | remaining {len(remaining)} | "
          f"{len(batches)} batches of <= {args.batch_size}", file=sys.stderr)

    result = {'batch_index': list(range(len(batches))),
              'n_remaining': len(remaining), 'n_batches': len(batches)}
    print(json.dumps(result))

    github_output = os.environ.get('GITHUB_OUTPUT')
    if github_output:
        with open(github_output, 'a') as f:
            f.write(f"matrix={json.dumps({'batch_index': result['batch_index']})}\n")
            f.write(f"count={len(batches)}\n")


if __name__ == '__main__':
    main()
