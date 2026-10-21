# tests

`fixtures/partials/v10/` holds the partial-sum parquets (partials schema 3) of two v10 tiles: 033_198 (Qin Ling: Huai He,
Yellow River, a mountain range, two water-tower parts) and 053_070 (Colombian Andes: Magdalena, San Juan), all five unit
types, built by `tests/build_fixtures.py` on 2026-10-20 (read-only: the ancillary and store windows, with the
`WTU_SIGNED` layer rasterized in memory before its fleet wave). They are the input of the smoke test in
`.github/workflows/ci.yml`: install the `ci` environment, import the package, and execute the three aggregation
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
the GTOPO30 histogram and the level-6 population table are tracked. Regenerate the fixtures with
`tests/build_fixtures.py <tiles>` after a schema change (or copy two tiles from `partials/<version>/` once the fleet has
re-mapped). To run the smoke test locally without your SAS token being picked up, put a clone of the production repo
without `config/sas_token.txt` first on `PYTHONPATH`.

`validate_schema3.py <old aggregate.py> <pixel tables...>` checks the schema-3 map and reduce against the schema-2 code
on real pixel tables (see `docs/aggregation_lineage.md`).
