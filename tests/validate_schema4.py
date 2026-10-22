"""Validate partials schema 4 (window-relative onset, aggregate.WINDOW_M).

1. Synthetic tile with a known answer: the west part of a 32 km tile melts 30 days earlier and is mostly
   north-facing, the east part mostly south-facing, and every pixel carries the same aspect effect
   (north-facing onset 10 days later than south-facing). The range-wide comparison
   (aggregate.elevation_relative, the schema-3 triplet panel) is biased by location. (a) With the boundary on a
   window edge the window-relative onset must recover N-S = +10 d, for the median and every water year.
   (b) With the boundary inside a window the window that straddles it keeps some of the bias; reported, not
   asserted (a sharp regime boundary is the worst case of the method).
2. Real tiles (needs the Azure token): the pixel table from datacube.tabulate_tile, mapped with the schema-3
   code (``--old <path to the schema-3 aggregate.py>``, e.g. from ``git show main:gsro_analysis/aggregate.py``)
   and with this code. Every schema-3 column must be bit-identical; reports map time and partials size, and
   the window-relative N-S / E-W per mountain range on the tile next to the range-wide value.

Usage: python tests/validate_schema4.py [--old OLD_AGGREGATE_PY] [ROW_COL ...]
"""
import argparse, importlib.util, os, sys, tempfile, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gsro_analysis import aggregate  # noqa: E402

WY = list(range(2015, 2026))


def synthetic_table(seed=0, size=400, true_ns=10.0, boundary_m=520_000):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:size, 0:size]
    x = 500_000 + xx.ravel() * 80.0 + 40; y = 4_000_000 + yy.ravel() * 80.0 + 40
    west = x < boundary_m
    north_share = np.where(west, 0.8, 0.2)                       # aspect make-up varies with location
    north = rng.random(x.size) < north_share
    aspect = np.where(north, rng.uniform(-60, 60, x.size) % 360, rng.uniform(120, 240, x.size))
    dem = rng.uniform(3000, 3199, x.size)                         # two 100 m bands
    base = np.where(west, 150.0, 180.0)                           # location effect
    med = base + true_ns / 2 * np.cos(np.deg2rad(aspect)) + rng.normal(0, 3, x.size)
    df = pd.DataFrame({"tile_row": 0, "tile_col": 0, "x": x.astype("float32"), "y": y.astype("float32"),
                       "original_lat": 36.0, "original_lon": 75.0,
                       "runoff_onset_median": np.round(med).astype("int16"), "temporal_resolution_median": 6.0,
                       "runoff_onset_mad": 5.0, "dem": dem.astype("int16"), "aspect": np.round(aspect).astype("int16") % 360,
                       "slope": 25, "chili": 0.5, "snow_classification": 2, "esa_worldcover": 60, "forest_cover_fraction": 0,
                       "GMBA_V2_ID": 1, "PFAF_ID": 4_000_000 + 1, "continent": 2, "GRDC_MRBID": -9999, "WTU_SIGNED": -9999})
    for j, wy in enumerate(WY):
        df[f"runoff_onset_WY{wy}"] = np.round(med + rng.normal(0, 4, x.size)).astype("int16")
    return df


def contrasts(ds, var, nvar):
    """pixel-weighted N-S / E-W over the unit's bands (mean over elevation of the aspect curve, weighted by n)"""
    w = ds[nvar].where(np.isfinite(ds[var]), 0)
    curve = (ds[var].fillna(0) * w).sum("elevation") / w.sum("elevation")
    c = aggregate.aspect_contrast(curve, w.sum("elevation"))
    return float(c.ns.squeeze()), float(c.ew.squeeze())


def test_synthetic(boundary_m=520_000, strict=True):
    df = synthetic_table(boundary_m=boundary_m)
    parts = aggregate.tile_partials(df, WY, unit_types=["mountain_ranges"])
    sparse = aggregate.sparse_from_frame(parts, "mountain_ranges")
    ds = aggregate.view(sparse, "mountain_ranges", ("elevation", "aspect"), WY)
    pooled = aggregate.elevation_relative(ds["runoff_onset_median"])
    ds["pooled"] = pooled
    ns_pooled, _ = contrasts(ds, "pooled", "runoff_onset_median_n")
    ns_w, ew_w = contrasts(ds, "runoff_onset_wrel", "runoff_onset_wrel_n")
    yearly = [contrasts(ds.sel(water_year=y), "runoff_onset_wrel_wy", "runoff_onset_wrel_wy_n")[0] for y in WY]
    share = float(ds.runoff_onset_wrel_n.sum() / ds.runoff_onset_median_n.sum())
    print(f"synthetic, boundary at {boundary_m / 1000:.0f} km ({'window edge' if boundary_m % aggregate.WINDOW_M == 0 else 'inside a window'}): true N-S +10.0 d | range-wide (schema-3 panel) {ns_pooled:+.2f} | window-relative {ns_w:+.2f} "
          f"(E-W {ew_w:+.2f}) | per year {min(yearly):+.2f} .. {max(yearly):+.2f} | pixels with a value {share:.0%}")
    if strict:
        assert abs(ns_w - 10) < 0.6, ns_w
        assert all(abs(v - 10) < 1.0 for v in yearly), yearly
    assert abs(ns_pooled - 10) > 5, "the synthetic case should bias the range-wide panel"
    return ns_pooled, ns_w


def test_forest_groups():
    """Half the pixels forested (FCF 80 %), melting 6 d earlier and mostly north-facing: each forest-cover group is
    referenced to itself, so the open (FCF <= 50 %) and forested (> 50 %) views both recover the +10 d aspect effect,
    and neither is biased by the other group's timing."""
    df = synthetic_table(boundary_m=520_000)
    rng = np.random.default_rng(3)
    forest = rng.random(len(df)) < np.where(np.cos(np.deg2rad(df.aspect.to_numpy())) > 0, 0.7, 0.3)
    df.loc[forest, "forest_cover_fraction"] = 80
    for c in ["runoff_onset_median"] + [f"runoff_onset_WY{y}" for y in WY]:
        df.loc[forest, c] = (df.loc[forest, c] - 6).astype("int16")
    parts = aggregate.tile_partials(df, WY, unit_types=["mountain_ranges"])
    sparse = aggregate.sparse_from_frame(parts, "mountain_ranges")
    res = {}
    for lab, fcf in (("open (FCF <= 50 %)", list(range(5))), ("forest (FCF > 50 %)", list(range(5, 10))), ("all", list(range(11)))):
        ds = aggregate.view(sparse, "mountain_ranges", ("elevation", "aspect"), WY, fcf=fcf)
        res[lab] = contrasts(ds, "runoff_onset_wrel", "runoff_onset_wrel_n")[0]
    print("forest groups: window-relative N-S " + " | ".join(f"{k} {v:+.2f}" for k, v in res.items()) + " (true +10.0)")
    assert all(abs(v - 10) < 0.8 for v in res.values()), res


def test_range_metrics():
    """range_metrics on the synthetic tile: the window method (default on a schema-4 cube) recovers the +10 d
    sunny-shady difference with the earliest onset on south-facing slopes, per range, per year and per band; the
    pooled method (schema 3) does not."""
    from gsro_analysis import range_metrics as rm
    df = synthetic_table(boundary_m=520_000)
    ds = aggregate.view(aggregate.sparse_from_frame(aggregate.tile_partials(df, WY, unit_types=["mountain_ranges"]), "mountain_ranges"),
                        "mountain_ranges", ("elevation", "aspect"), WY)
    ds = ds.assign_coords(mountain_range=["synthetic"])
    v = rm.base_view(ds)
    clim_w, yr_w, _ = rm.sunny_shady(v); clim_p, _, _ = rm.sunny_shady(v, method="pooled")
    lapse, amp = rm.yearly_topography(v)
    band = rm.band_aspect(v)
    bm = band[band.water_year == "median"]
    print(f"range_metrics: sunny-shady window {clim_w.sunny_shady_d.iloc[0]:.2f} d at {clim_w.earliest_aspect_deg.iloc[0]:.0f} deg "
          f"(pooled {clim_p.sunny_shady_d.iloc[0]:.2f} d at {clim_p.earliest_aspect_deg.iloc[0]:.0f} deg) | per year "
          f"{yr_w.min(axis=1).iloc[0]:.2f}..{yr_w.max(axis=1).iloc[0]:.2f} | yearly_topography amp {amp.mean(axis=1).iloc[0]:.2f} | "
          f"bands 2A {2 * bm.amp_d.min():.2f}..{2 * bm.amp_d.max():.2f}")
    assert abs(clim_w.sunny_shady_d.iloc[0] - 10) < 0.8 and abs(clim_w.earliest_aspect_deg.iloc[0] - 180) < 15
    assert (yr_w.sub(10).abs() < 1.2).all(axis=None) and (amp.sub(10).abs() < 1.2).all(axis=None)
    assert (bm.amp_d.mul(2).sub(10).abs() < 1.5).all()


def load_module(path, name="aggregate_old"):
    spec = importlib.util.spec_from_file_location(name, path); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def real_tiles(tiles, old_path):
    try: open('/proc/self/oom_score_adj', 'w').write('1000')
    except OSError: pass
    from gsro_analysis import datacube, settings
    old = load_module(old_path)
    config = settings.load_config(); repo = datacube.open_ancillary_repo(config)
    gds = config.open_runoff_onset_dataset(chunks=None, mask_and_scale=True)
    m = pd.read_csv(Path(__file__).resolve().parents[1] / "analyses/mountain_ranges/results/v10/mountain_range_metrics.csv").set_index("GMBA_V2_ID")
    rows = []
    for t in tiles:
        r, c = map(int, t.split("_"))
        t0 = time.time(); df = datacube.tabulate_tile(config, r, c, datacube.open_ancillary_window(config, r, c, repo), global_ds=gds); t_tab = time.time() - t0
        t0 = time.time(); p_old = old.tile_partials(df, config.water_years); t_old = time.time() - t0
        t0 = time.time(); p_new = aggregate.tile_partials(df, config.water_years); t_new = time.time() - t0
        keys = [k for k in aggregate.KEY_COLS]
        a = p_old.sort_values(keys, na_position="first").reset_index(drop=True)
        b = p_new.sort_values(keys, na_position="first").reset_index(drop=True)
        same_rows = len(a) == len(b)
        common = [col for col in a.columns if col in b.columns]
        identical = same_rows and all(a[col].equals(b[col]) or np.array_equal(a[col].to_numpy(), b[col].to_numpy(), equal_nan=True) for col in common)
        new_cols = [col for col in b.columns if col not in a.columns]
        with tempfile.TemporaryDirectory() as d:
            fo, fn = Path(d) / "o.parquet", Path(d) / "n.parquet"
            p_old.to_parquet(fo, index=False); p_new.to_parquet(fn, index=False); so, sn = os.path.getsize(fo), os.path.getsize(fn)
        print(f"tile {t}: {len(df):,} px | tabulate {t_tab:.0f}s, map old {t_old:.1f}s new {t_new:.1f}s | rows {len(a)} -> {len(b)} | "
              f"schema-3 columns identical: {identical} | +{len(new_cols)} columns | parquet {so / 1e6:.2f} -> {sn / 1e6:.2f} MB", flush=True)
        assert identical, f"tile {t}: schema-3 columns changed"
        w0 = aggregate.WINDOW_M; aggregate.WINDOW_M = 5_000                    # the same at 5 km windows, for scale
        try: p5 = aggregate.tile_partials(df, config.water_years, unit_types=["mountain_ranges"])
        finally: aggregate.WINDOW_M = w0
        ds5 = aggregate.view(aggregate.sparse_from_frame(p5, "mountain_ranges"), "mountain_ranges", ("elevation", "aspect"), config.water_years)
        del df
        sparse = aggregate.sparse_from_frame(p_new, "mountain_ranges")
        ds = aggregate.view(sparse, "mountain_ranges", ("elevation", "aspect"), config.water_years)
        ds["pooled"] = aggregate.elevation_relative(ds["runoff_onset_median"])
        for uid in ds.mountain_range.values:
            d = ds.sel(mountain_range=uid)
            if int(d.runoff_onset_wrel_n.sum()) < 20_000:
                continue
            ns_w, ew_w = contrasts(d, "runoff_onset_wrel", "runoff_onset_wrel_n")
            ns_p, ew_p = contrasts(d, "pooled", "runoff_onset_median_n")
            ns_5, _ = contrasts(ds5.sel(mountain_range=uid), "runoff_onset_wrel", "runoff_onset_wrel_n") if uid in ds5.mountain_range else (np.nan, np.nan)
            rows.append({"tile": t, "GMBA_V2_ID": int(uid), "name": m["name"].get(int(uid), str(uid)), "px_wrel": int(d.runoff_onset_wrel_n.sum()),
                         "share_with_wrel": float(d.runoff_onset_wrel_n.sum() / d.runoff_onset_median_n.sum()),
                         "NS_window": ns_w, "NS_window_5km": ns_5, "EW_window": ew_w, "NS_tile_pooled": ns_p, "EW_tile_pooled": ew_p})
    out = pd.DataFrame(rows); print(out.round(2).to_string(index=False)); return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--old"); ap.add_argument("--out"); ap.add_argument("tiles", nargs="*")
    args = ap.parse_args()
    test_synthetic(520_000, strict=True)
    test_synthetic(516_000, strict=False)
    test_forest_groups()
    test_range_metrics()
    if args.tiles:
        res = real_tiles(args.tiles, args.old)
        if args.out: res.to_csv(args.out, index=False)
