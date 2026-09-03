# tests

`fixtures/aggregated_results/v10/partials/` holds the partial-sum parquets of two v10 dry-run tiles
(001_183: basins and a continent only; 094_077: two mountain ranges as well), written by the fleet map
on 2026-09-03 with the current schema (HydroBASINS level-6 basin ids, the continents aspect key, the
`fcf_lte_50` and `full_dataset` filter tags). They are the input of the smoke test in
`.github/workflows/ci.yml`: install the `ci` environment, import the package, and run the reduce on
them without any credential:

```bash
GSRO_AGGREGATED_ROOT=tests/fixtures/aggregated_results \
  pixi run -e ci python pipeline/scripts/reduce_partials.py --no-download
```

which must write the three default cubes under `tests/fixtures/aggregated_results/v10/`. The GMBA
and continents polygons are downloaded into `data/geometries/sources/` on first use (≈ 34 MB); the
GTOPO30 land histogram is tracked (`data/gtopo30_lat_elev_histogram.nc`). Regenerate the fixtures by
copying two tiles from the local partials cache after a schema change.
