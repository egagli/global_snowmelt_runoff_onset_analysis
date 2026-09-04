# tests

`fixtures/partials/v10/` holds the partial-sum parquets of two v10 dry-run tiles (001_183: basins and a continent only;
094_077: two mountain ranges as well), written by the fleet map on 2026-09-03 with the current schema (HydroBASINS level-6
basin ids, the continents aspect key, the `fcf_lte_50` and `full_dataset` filter tags). They are the input of the smoke
test in `.github/workflows/ci.yml`: install the `ci` environment, import the package, and execute the three aggregation
notebooks headlessly on them without any credential (the partials download and the ERA5-Land step skip themselves when
no SAS token is present):

```bash
GSRO_PARTIALS_ROOT=tests/fixtures/partials GSRO_OUTPUT_ROOT=/tmp/ci/analyses \
  pixi run -e ci jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=python3 \
    --output-dir /tmp/ci/executed analyses/continents/0_aggregate_by_continent.ipynb \
    analyses/mountain_ranges/0_aggregate_by_mountain_range.ipynb analyses/river_basins/0_aggregate_by_river_basin.ipynb
```

which must write the cubes under `/tmp/ci/analyses/<unit>/data/aggregation/v10/` and the two metrics tables under
`/tmp/ci/analyses/<unit>/results/v10/`. GMBA and the USGS continents are read from the web; the BasinATLAS gdb (2.7 GB)
is downloaded into `analyses/river_basins/data/geometries/` on first use and restored from the Actions cache afterwards;
the GTOPO30 histogram and the level-6 population table are tracked. Regenerate the fixtures by copying two tiles from
`partials/<version>/` after a schema change.
