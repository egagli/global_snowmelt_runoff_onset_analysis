"""Write the CI fixture partials (the current PARTIALS_SCHEMA) for a few tiles, read-only: the ancillary window from the live store,
with any unit-id layer the store does not hold yet (WTU_SIGNED before its fleet wave) rasterized in memory exactly
as the ancillary job would; the store window; then tabulate_tile -> tile_partials -> tests/fixtures/partials/<version>/.

    GSRO env with the SAS token: pixi run -e ci python tests/build_fixtures.py 033_198 053_070"""
import sys
import numpy as np
from pathlib import Path
from gsro_analysis import settings, datacube, aggregate

cfg = settings.load_config(); repo = datacube.open_ancillary_repo(cfg)
out = Path(__file__).parent / 'fixtures' / 'partials' / cfg.version; out.mkdir(parents=True, exist_ok=True)
for rc in sys.argv[1:]:
    r, c = map(int, rc.split('_'))
    win = datacube.open_ancillary_window(cfg, r, c, repo)
    missing = [l for l in datacube.UNIT_ID_LAYERS if l not in win]
    if missing:
        ds = datacube.add_mountain_range_and_basin_and_continent(cfg.get_tile(r, c), win[['dem']].copy(), layers=missing)
        for l in missing:
            v = ds[l].astype('float32')
            win[l] = v.where(v != -9999)
    df = datacube.tabulate_tile(cfg, r, c, win)
    parts = aggregate.tile_partials(df, cfg.water_years)
    parts.insert(0, 'tile_col', np.int16(c)); parts.insert(0, 'tile_row', np.int16(r))
    path = out / f'tile_{r:03d}_{c:03d}.parquet'
    datacube.write_partials_parquet(parts, path)
    print(rc, f'{len(df):,} pixels, {len(parts):,} rows ({sorted(parts.unit_type.unique())}), added in memory: {missing} -> {path}', flush=True)
