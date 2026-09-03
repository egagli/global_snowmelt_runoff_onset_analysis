"""The ERA5-Land store operations — the ERA5 Acquire workflow's entry point, equally
runnable locally (see gsro_analysis/era5.py for the layout and the ledger).

  era5_land.py init    [--start-fresh]      create the version's repository with the empty template
                                             if it does not exist (--start-fresh: DELETE and recreate)
  era5_land.py plan    [--manifest-out F]   fold the commit ledger; print the remaining water years as
                                             json (with $GITHUB_OUTPUT: matrix=, count=, anomaly_needed=)
  era5_land.py fetch   --water-year 2020    fetch one water year from Earth Engine and commit it
  era5_land.py anomaly [--force]            build the anomaly group once every water year is committed
                                             (a no-op while years are missing or when it is up to date)
  era5_land.py status                        print the ledger

--config picks the dataset version (default settings.CONFIG_FILE); --local-store PATH
runs against a local icechunk repository instead of Azure (tests). Every water year of
config.water_years is acquired from Earth Engine; nothing is ever copied from another
version. Needs the Azure SAS token; fetch needs Earth Engine.
"""

import argparse
import json
import os
import sys
import time

from gsro_analysis import era5, settings


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('command', choices=['init', 'plan', 'fetch', 'anomaly', 'status'])
    p.add_argument('--config', default=settings.CONFIG_FILE)
    p.add_argument('--local-store', help='path of a local icechunk repository (tests)')
    p.add_argument('--start-fresh', action='store_true',
                   help='init: delete an existing repository and recreate the empty template')
    p.add_argument('--water-year', type=int, help='fetch: the water year to acquire')
    p.add_argument('--force', action='store_true', help='anomaly: rebuild even if up to date')
    p.add_argument('--manifest-out', default='era5_plan.json', help='plan: where to write the json')
    p.add_argument('--ee-key', help='Earth Engine service-account key json '
                   '(default: the production clone\'s config/ee_key.json)')
    return p.parse_args()


def print_status(config, st):
    years = [int(y) for y in config.water_years]
    done = sorted(st['years'])
    print(f"{era5.repo_prefix(config)} @ {st['snapshot_id']}")
    print(f"water years {years[0]}-{years[-1]}: {len(done)} committed {done} | remaining {st['remaining']}")
    if st['anomaly'] is None:
        print("anomaly: none")
    else:
        a = st['anomaly']
        print(f"anomaly: commit {a['snapshot_id']} ({a['written_at']}), base {a.get('base_water_years')}"
              f"{' — STALE' if st['anomaly_stale'] else ''}")


def main():
    args = parse_args()
    config = settings.load_config(args.config)
    log = lambda m: print(m, flush=True)  # noqa: E731

    if args.command == 'init':
        era5.initialize(config, start_fresh=args.start_fresh, local_store=args.local_store, log=log)
        return

    repo = era5.open_repo(config, args.local_store)

    if args.command == 'status':
        print_status(config, era5.status(config, repo))

    elif args.command == 'plan':
        st = era5.status(config, repo)
        print_status(config, st)
        anomaly_needed = st['anomaly'] is None or st['anomaly_stale']
        plan = {'config': args.config, 'version': config.version, 'snapshot_id': st['snapshot_id'],
                'water_years': [int(y) for y in config.water_years], 'remaining': st['remaining'],
                'anomaly_needed': anomaly_needed}
        with open(args.manifest_out, 'w') as f:
            json.dump(plan, f)
        print(json.dumps({'remaining': st['remaining'], 'anomaly_needed': anomaly_needed}))
        github_output = os.environ.get('GITHUB_OUTPUT')
        if github_output:
            with open(github_output, 'a') as f:
                f.write(f"matrix={json.dumps({'year': st['remaining']})}\n")
                f.write(f"count={len(st['remaining'])}\n")
                f.write(f"anomaly_needed={'true' if anomaly_needed else 'false'}\n")

    elif args.command == 'fetch':
        if args.water_year is None:
            sys.exit('fetch requires --water-year')
        if args.water_year not in [int(y) for y in config.water_years]:
            sys.exit(f"WY{args.water_year} is not in config.water_years")
        print(f"Earth Engine initialized as {settings.initialize_earthengine(key_file=args.ee_key)}", flush=True)
        t0 = time.time()
        era5.acquire_water_year(config, args.water_year, repo=repo, log=log)
        log(f"WY{args.water_year} done in {time.time() - t0:.0f}s")

    elif args.command == 'anomaly':
        era5.build_anomaly(config, repo=repo, force=args.force, log=log)


if __name__ == '__main__':
    main()
