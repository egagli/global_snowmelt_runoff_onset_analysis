"""Population per HydroBASINS level-6 basin -> the tracked table
analyses/river_basins/data/geometries/hydrobasins_level6_population.csv (HYBAS_ID, PFAF_ID, total_population).

Source: GPW v4.11 population count (CIESIN, the most recent epoch, native 30 arc-second grid,
settings.GPW_POPULATION_COLLECTION) summed on Earth Engine inside every polygon of HydroBASINS v1
level 6 (settings.HYDROBASINS_LEVEL6_ASSET; the level the fleet keys the basin partials by). This is
the Python-API port of the April 2025 Code Editor export that produced the level-5 geojson; level 6 is
finer, and because Pfafstetter codes nest by digit prefix the level-5 population is the exact sum of
the level-6 rows (PFAF_ID // 10), which is what analyses/river_basins/0_aggregate_by_river_basin.ipynb
does. Geometry is not exported: the notebooks take the polygons from the BasinATLAS gdb they already
read, so the table is ~0.5 MB and lives in git.

Runs in chunks of --chunk-size basins (reduceRegions with tileScale 16, geometry dropped from the
result), writes every finished chunk to <out>.part.csv and resumes from it, so an interrupted run
loses at most one chunk. Needs the Earth Engine service key (settings.initialize_earthengine()).
~16,400 basins take on the order of 15-30 minutes.

    pixi run population                       # = python pipeline/scripts/get_basin_population.py
"""

import argparse
import csv
import datetime as dt
import time
from pathlib import Path

from gsro_analysis import paths, settings

OUT_NAME = 'hydrobasins_level6_population.csv'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=str(paths.geometries('river_basins') / OUT_NAME))
    p.add_argument('--chunk-size', type=int, default=200, help='basins per reduceRegions request')
    p.add_argument('--force', action='store_true', help='rebuild even if the output exists')
    p.add_argument('--ee-key', help='Earth Engine service-account key json (default: the production clone\'s config/ee_key.json)')
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.out)
    if out.exists() and not args.force:
        print(f"{out} exists; --force rebuilds it")
        return
    import ee
    print(f"Earth Engine initialized as {settings.initialize_earthengine(key_file=args.ee_key)}", flush=True)

    population_image = (ee.ImageCollection(settings.GPW_POPULATION_COLLECTION)
                        .sort('system:time_start', False).first().select('population_count'))
    epoch = ee.Date(population_image.get('system:time_start')).format('YYYY').getInfo()
    scale_m = population_image.projection().nominalScale().getInfo()
    basins = ee.FeatureCollection(settings.HYDROBASINS_LEVEL6_ASSET)
    hybas_ids = basins.aggregate_array('HYBAS_ID').getInfo()
    print(f"GPW v4.11 epoch {epoch} at {scale_m:.0f} m | {len(hybas_ids)} level-6 basins", flush=True)

    # resume: rows already written to the .part file
    part = out.with_name(out.stem + '.part.csv')
    done = {}
    if part.exists():
        with open(part) as f:
            for row in csv.DictReader(f):
                done[int(row['HYBAS_ID'])] = row
        print(f"resuming: {len(done)} basins already in {part.name}", flush=True)
    todo = [h for h in hybas_ids if h not in done]
    chunks = [todo[i:i + args.chunk_size] for i in range(0, len(todo), args.chunk_size)]

    t0 = time.time()
    with open(part, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['HYBAS_ID', 'PFAF_ID', 'total_population'])
        if not done:
            writer.writeheader()
        for k, chunk in enumerate(chunks):
            fc = basins.filter(ee.Filter.inList('HYBAS_ID', chunk))
            summed = population_image.reduceRegions(collection=fc, reducer=ee.Reducer.sum(),
                                                    scale=scale_m, tileScale=16)
            summed = summed.select(['HYBAS_ID', 'PFAF_ID', 'sum'], None, False)   # drop the geometry
            for attempt in range(1, 6):
                try:
                    features = summed.getInfo()['features']
                    break
                except Exception as e:  # noqa: BLE001 - EE quota / transient
                    if attempt == 5:
                        raise
                    delay = 30 * attempt
                    print(f"  chunk {k}: {type(e).__name__}: {e}; retry {attempt}/5 in {delay}s", flush=True)
                    time.sleep(delay)
            for feat in features:
                props = feat['properties']
                writer.writerow({'HYBAS_ID': int(props['HYBAS_ID']), 'PFAF_ID': int(props['PFAF_ID']),
                                 'total_population': round(float(props.get('sum') or 0.0))})
            f.flush()
            print(f"  chunk {k + 1}/{len(chunks)} ({len(chunk)} basins) done ({time.time() - t0:.0f}s)", flush=True)

    # the final table: sorted, with two provenance comment lines (read with pd.read_csv(..., comment='#'))
    with open(part) as f:
        rows = sorted(csv.DictReader(f), key=lambda r: int(r['PFAF_ID']))
    total = sum(float(r['total_population']) for r in rows)
    with open(out, 'w', newline='') as f:
        f.write(f"# HydroBASINS v1 level-6 population: {settings.GPW_POPULATION_COLLECTION} (epoch {epoch}, "
                f"{scale_m:.0f} m) summed per {settings.HYDROBASINS_LEVEL6_ASSET} polygon on Earth Engine\n")
        f.write(f"# generated {dt.date.today().isoformat()} by pipeline/scripts/get_basin_population.py; "
                f"level 5 = groupby(PFAF_ID // 10).sum()\n")
        writer = csv.DictWriter(f, fieldnames=['HYBAS_ID', 'PFAF_ID', 'total_population'])
        writer.writeheader()
        writer.writerows(rows)
    part.unlink()
    print(f"wrote {out}: {len(rows)} basins, total population {total / 1e9:.3f} billion ({time.time() - t0:.0f}s)")


if __name__ == '__main__':
    main()
