"""Schema-3 validation: on real pixel tables, the schema-3 map + sparse reduce + cube() must reproduce the
schema-2 cubes (from the schema-2 code at the parent commit) once their chili_class dimension is collapsed:
pixel counts exactly, means and stds to float32 precision, for every default group and both former filter tags.
Also checks the DuckDB reduce against the in-memory reduce, the FCF classes against the old FCF filter, and that
terrain counts bound the analysed counts.

    pixi run -e ci python tests/validate_schema3.py <old_aggregate.py> <pixel_table.parquet> [...]

Pixel tables must hold every pixel with a DEM (schema-3 tabulate_tile) and a WTU_SIGNED column."""
import sys, importlib.util, tempfile, numpy as np, pandas as pd, xarray as xr
from pathlib import Path
from gsro_analysis import aggregate as new, datacube

spec = importlib.util.spec_from_file_location("old_aggregate", sys.argv[1]); old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
WY = list(range(2015, 2026))
pix = [pd.read_parquet(f) for f in sys.argv[2:]]
onset_only = [p[p.runoff_onset_median > 0] for p in pix]
old_parts = pd.concat([old.tile_partials(p, WY) for p in onset_only], ignore_index=True)
new_parts = pd.concat([new.tile_partials(p, WY) for p in pix], ignore_index=True)
print(f"rows: schema 2 {len(old_parts):,} | schema 3 {len(new_parts):,} ({int(new_parts.onset_bin.notna().sum()):,} histogram)")

COUNTS = ["runoff_onset_median_n", "runoff_onset_mad_n", "runoff_onset_n", "temporal_resolution_median_n"]
MEANS = ["runoff_onset_median", "runoff_onset_median_std", "runoff_onset_mad", "runoff_onset", "runoff_onset_anomaly",
         "runoff_onset_anomaly_std", "temporal_resolution_median", "n_years", "fcf_corr"]
fail = 0
# continents_aspect is left out: its dense lat x elev x aspect x CHILI x year reference cube needs > 8 GB to align
for group in ("mountain_ranges", "river_basins", "river_basins_l6", "continents", "grdc_basins", "grdc_basins_aspect"):
    for tag in ("fcf_lte_50", "full_dataset"):
        try:
            a = old.collapse(old.reduce_partials(old_parts, group, tag, WY))
        except ValueError as e:
            print(f"  {group}/{tag}: schema 2 has no rows ({e})"); continue
        b = new.reduce_partials(new_parts, group, tag, WY)
        a, b = xr.align(a, b, join="outer")
        msgs = []
        for v in COUNTS:
            d = int(np.abs(a[v].fillna(0) - b[v].fillna(0)).max())
            if d: msgs.append(f"{v} max |diff| {d}")
        for v in MEANS:
            x, y = a[v].values, b[v].values
            both = np.isfinite(x) & np.isfinite(y)
            if (np.isfinite(x) != np.isfinite(y)).any(): msgs.append(f"{v} NaN pattern differs ({int((np.isfinite(x) != np.isfinite(y)).sum())})")
            if both.any():
                # stds: the schema-2 reference pools float32 per-CHILI stds (cancellation in ss/n - mean^2), so an
                # absolute tolerance of 0.01 d; means: relative 2e-5
                if v.endswith("_std"):
                    dd = np.nanmax(np.abs(x[both] - y[both]))
                    if dd > 0.01: msgs.append(f"{v} max |diff| {dd:.3f} d")
                else:
                    rel = np.nanmax(np.abs(x[both] - y[both]) / np.maximum(1, np.abs(x[both])))
                    if rel > 2e-5: msgs.append(f"{v} max rel diff {rel:.2e}")
        fail += bool(msgs)
        print(f"  {group}/{tag}: {'OK' if not msgs else 'DIFF ' + '; '.join(msgs)}  (cells {int((b.runoff_onset_median_n > 0).sum()):,})")

# DuckDB reduce == in-memory reduce
with tempfile.TemporaryDirectory() as d:
    paths = []
    for i, p in enumerate(pix):
        f = Path(d) / f"t{i}.parquet"; datacube.write_partials_parquet(new.tile_partials(p, WY), f); paths.append(f)
    for ut in new.UNIT_TYPES:
        for h in (False, True):
            if h and ut not in new.HIST_UNIT_TYPES: continue
            s1 = new.reduce_sparse(paths, ut, histograms=h)
            s2 = new.sparse_from_frame(new_parts, ut, histograms=h)
            k = [c for c in new.KEY_COLS if c in s1]
            s1 = s1.sort_values(k, na_position="last").reset_index(drop=True); s2 = s2.sort_values(k, na_position="last").reset_index(drop=True)
            cols = [c for c in s2.columns if c not in new.KEY_COLS]
            ok = len(s1) == len(s2) and np.allclose(s1[cols].fillna(0).to_numpy(float), s2[cols].fillna(0).to_numpy(float), rtol=1e-12, atol=1e-9)
            fail += not ok
            print(f"  duckdb vs pandas {ut}{' hist' if h else ''}: {'OK' if ok else 'DIFF'} ({len(s1):,} vs {len(s2):,} rows)")
# terrain bounds analysed; histograms sum to the FCF<=50 analysed pixel counts
s = pd.concat([new.sparse_from_frame(new_parts, ut) for ut in new.UNIT_TYPES])
bad = int((s.n > s.n_terrain).sum()); fail += bad > 0
print(f"  rows with analysed n > n_terrain: {bad}")
for ut in new.HIST_UNIT_TYPES:
    h = new.sparse_from_frame(new_parts, ut, histograms=True); t = new.sparse_from_frame(new_parts, ut)
    m = int(t[t.fcf_class <= 4].n.sum()); hm = int(h.h_med.sum()); fail += m != hm
    print(f"  {ut}: histogram median counts {hm:,} vs analysed FCF<=50 pixels {m:,} {'OK' if m == hm else 'DIFF'}")
print("FAILURES:", fail); sys.exit(1 if fail else 0)
