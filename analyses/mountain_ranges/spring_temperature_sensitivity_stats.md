# Spring-temperature sensitivity of runoff onset: the statistics, what they mean, and how to report them

Status: written 2026-09-08 against dataset version **v10** (water years 2015–2025, 11 years). Every number
below was produced by the code blocks in this document, run in that order on the v10 mountain-range cube
(`data/aggregation/v10/all_mountain_ranges_fcf_lte_50.nc`) and the v10 metrics table
(`results/v10/mountain_range_metrics.csv`). The per-range statistics (sections 1–3) are what
`0_aggregate_by_mountain_range.ipynb` writes into the metrics table and what the two composite world maps
draw; the aggregate tests (sections 4–8) are **not yet in a notebook**. They live in this document and in a
scratch script, and should be ported into `temperature_sensitivity.ipynb` before anything is quoted from
them in a publication.

The question every section serves: *do warmer springs bring earlier snowmelt runoff onset, how much
earlier per degree, and how confident can we be with at most eleven years per mountain range?*

---

## Glossary

Terms in the order they are needed. Everything else in this document is plain arithmetic on these.

| term | meaning here |
| --- | --- |
| **water year (WY)** | 1 Oct–30 Sep in the northern hemisphere, 1 Apr–31 Mar in the southern; WY2015 is the first, WY2025 the last (11 years) |
| **runoff onset** | the day of the water year (DOWY) on which snowmelt runoff begins at a pixel, from the Sentinel-1 backscatter minimum; the dataset's variable |
| **onset anomaly** | a pixel's onset in one water year minus its median onset over all water years, in days; negative = earlier than usual |
| **range mean anomaly** | the pixel-count-weighted mean of the onset anomaly over all valid (elevation × aspect × insolation-class) bins of a GMBA mountain range in one water year; one number per range per year |
| **spring temperature anomaly** | the range's mean ERA5-Land 2 m air temperature anomaly over the three spring months (MAM north, SON south), in °C, relative to the per-pixel median over the same water years |
| **GMBA** | Global Mountain Biodiversity Assessment mountain inventory v2.0, "standard 300" polygons: the mountain-range units |
| **ERA5-Land** | ECMWF's land-surface reanalysis at 0.1°; the source of the temperature anomalies |
| **n** | for one range, the number of water years with both anomalies (3 to 11) |
| **OLS** | ordinary least squares: the straight line minimising the squared vertical misfits; its slope is "days of onset anomaly per °C of spring anomaly" |
| **Pearson r** | linear correlation coefficient, −1 to 1; r² is the share of variance the line explains |
| **p-value** | the probability of a statistic at least this extreme if the null hypothesis (no relationship) were true; p < 0.05 is the conventional "significant" |
| **Spearman ρ (rho)** | Pearson r computed on ranks: measures monotonic association, insensitive to one extreme point |
| **Kendall τ (tau)** | another rank correlation: the excess of concordant over discordant pairs; identical in spirit to ρ |
| **Theil–Sen slope** | the median of the slopes of all pairs of points (55 pairs for n = 11); a robust alternative to the OLS slope, with a distribution-free 95 % confidence interval (CI) from the same pairwise slopes |
| **CI** | confidence interval: the range of values compatible with the data at the stated confidence (95 % throughout) |
| **FDR, BH** | false discovery rate, controlled with the Benjamini–Hochberg procedure: when many tests are run, the share of "significant" results that are false; q < 0.05 means the expected share of false discoveries among those flagged is below 5 % |
| **IQR** | interquartile range, the 25th to 75th percentile |
| **pooled (fixed-effects) regression** | one regression through all range-years at once; because both anomalies are already centred within each range, this is equivalent to a regression with a separate intercept per range, i.e. it measures the *within-range* relationship |
| **block bootstrap (by year)** | resampling whole water years with replacement (all 160 ranges' values for that year move together) and refitting; the spread of the refits is the uncertainty that respects the fact that ranges share years |
| **permutation test** | rebuilding the statistic many times after shuffling the labels so the null hypothesis is true by construction; the share of shuffles that beat the observed value is the p-value |
| **field significance** | whether a *collection* of local tests (one per range) is collectively more significant than chance, given that the locations are not independent |
| **leave-one-year-out** | recomputing everything with one water year removed, once per year; shows whether a single year drives the result |
| **MAD** | median absolute deviation of onset across years (the dataset's interannual-variability statistic); not used below but appears in the same metrics table |

---

## 0. The data every statistic is computed from

**How.** For every range and water year, the range mean onset anomaly and the range mean spring temperature
anomaly, under the same masking rules the metrics table uses: bins with fewer than 100 pixels are
dropped; a range-year counts only if at least 10 % of the range's median-pixel count has data that year; a
range enters the analysis only with at least three paired years. That gives 160 ranges and 1580 range-years
(106 ranges have all 11 years; 20 have 7 or fewer).

**Why these rules.** Thin bins and thin years produce anomalies dominated by a handful of pixels; the
thresholds are the metrics table's, so every number here is on the same footing as the maps.

```python
import numpy as np, pandas as pd, xarray as xr
from scipy import stats as sps
from gsro_analysis import aggregate, paths

VERSION = 'v10'
cube = xr.open_dataset(paths.aggregation_dir('mountain_ranges', VERSION) / 'all_mountain_ranges_fcf_lte_50.nc')
BIN_DIMS = ['elevation', 'aspect', 'chili_class']
ds = aggregate.threshold(cube, 100)                                            # bins with < 100 pixels -> NaN
onset_anomaly = aggregate.weighted_mean(ds, 'runoff_onset_anomaly', BIN_DIMS)  # pixel-weighted range mean per water year
valid_fraction = ds['runoff_onset_n'].sum(BIN_DIMS) / ds['runoff_onset_median_n'].sum(BIN_DIMS)
onset_anomaly = onset_anomaly.where(valid_fraction >= 0.1)                     # a range-year needs 10 % of the range's pixels
spring_t = ds['temperature_2m'].sel(month=['spring_month_1', 'spring_month_2', 'spring_month_3']).mean('month')

pairs = pd.DataFrame({
    'range': np.repeat(onset_anomaly['mountain_range'].values, onset_anomaly.sizes['water_year']),
    'water_year': np.tile(onset_anomaly['water_year'].values, onset_anomaly.sizes['mountain_range']),
    'x': spring_t.transpose('mountain_range', 'water_year').values.ravel(),    # spring 2 m temperature anomaly [degC]
    'y': onset_anomaly.values.ravel(),                                         # range mean onset anomaly [days]
    'lat': np.repeat(ds['centroid_latitude'].values, onset_anomaly.sizes['water_year'])}).dropna()
n_pairs = pairs.groupby('range').size()
pairs = pairs[pairs['range'].isin(n_pairs[n_pairs >= 3].index)]               # a fit needs >= 3 paired years
years = np.sort(pairs['water_year'].unique())
```

| years per range (n) | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ranges | 1 | 7 | 8 | 10 | 5 | 13 | 10 | 106 |

Across all range-years the spring temperature anomaly has a standard deviation of 1.0 °C and the onset
anomaly 11.4 days.

---

## 1. Per-range OLS slope, Pearson r and p (what the main map shows)

**How.** For each range, OLS of the onset anomaly on the spring temperature anomaly over its n years:
slope [days/°C], r, p. This is the `anomaly_slope`, `anomaly_corr`, `anomaly_pval`, `anomaly_n` block of
the metrics table; the choropleth colour is the slope, the transparency is r.

**Why.** The slope is the quantity of interest in physical units (days earlier per degree warmer), r says
how much of the year-to-year scatter the line accounts for, and OLS is the estimator every reader knows.

```python
def per_range_fits(df):
    rows = []
    for name, g in df.groupby('range'):
        if len(g) < 3:
            continue
        ols = sps.linregress(g['x'], g['y'])
        ts = sps.mstats.theilslopes(g['y'], g['x'], 0.95)                     # section 3
        rho = sps.spearmanr(g['x'], g['y'])                                    # section 3
        rows.append(dict(range=name, n=len(g), lat=g['lat'].iloc[0], slope=ols.slope, r=ols.rvalue, p=ols.pvalue,
                         ts_slope=ts[0], ts_low=ts[2], ts_high=ts[3], rho=rho.statistic, rho_p=rho.pvalue))
    return pd.DataFrame(rows)

fits = per_range_fits(pairs)
```

**What it says (v10).**

| statistic over the 160 ranges | value |
| --- | --- |
| share with a negative slope (warmer → earlier) | 0.86 (137 of 160) |
| OLS slope: 5th / 25th / median / 75th / 95th percentile | −12.5 / −7.4 / **−4.9** / −2.0 / +3.8 days/°C |
| Pearson r: 5th / 25th / median / 75th / 95th percentile | −0.89 / −0.75 / **−0.60** / −0.18 / +0.37 |
| share with \|r\| > 0.4 | 0.68 (0.64 with r < −0.4) |
| ranges with a significant *positive* slope (p < 0.05) | 1 (Verkhoyansk Range, +3.1 days/°C, r = 0.64) |

By n, the relationship sharpens as the record fills in, which is what a real signal does and noise does not:

| n | ranges | median slope | median r | share p < 0.05 |
| --- | --- | --- | --- | --- |
| 5 | 7 | −0.5 | −0.07 | 0.00 |
| 6 | 8 | −3.5 | −0.58 | 0.38 |
| 7 | 10 | −3.5 | −0.54 | 0.10 |
| 8 | 5 | −1.9 | −0.12 | 0.00 |
| 9 | 13 | −2.0 | −0.26 | 0.23 |
| 10 | 10 | −7.2 | −0.50 | 0.20 |
| 11 | 106 | −5.8 | −0.67 | 0.61 |

**Significance for the study.** Per-range numbers are the *map*: the geography of sensitivity, which ranges
respond strongly (mid-latitude western North America, the Qilian Shan, the Dry Andes at −8 to −13 days/°C)
and which barely respond (Brooks Range, Kamchatka, the tropical Andes). They are descriptive. With n = 11 the
sampling uncertainty on one range's r is about ±0.2 and on its slope about ±7 days/°C (section 3), so a single
range's value is a rough estimate, not a measurement.

**Why report it.** It is the figure. Report the distribution across ranges (median, IQR, 5–95 % range) in the
text, and the per-range values with n and r in a supplementary table. Do **not** build a headline on a
threshold like "|r| > 0.4": the critical |r| for p < 0.05 at n = 11 is 0.60 (0.63 at n = 10, 0.71 at n = 8,
0.81 at n = 6), so |r| > 0.4 is below what eleven points can distinguish from noise.

---

## 2. Per-range significance and the multiple-testing correction

**How.** Count ranges with OLS p < 0.05, then adjust the 160 p-values with Benjamini–Hochberg and count
q < 0.05.

**Why.** 160 tests at the 5 % level would flag 8 ranges by chance alone; the FDR adjustment says how many
survive that accounting.

```python
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.minimum(q, 1)
    return out

fits['q'] = bh_fdr(fits['p'])
```

| | count of 160 |
| --- | --- |
| OLS p < 0.05 | 74 |
| OLS p < 0.10 | 87 |
| BH q < 0.05 | 47 |
| BH q < 0.10 | 73 |
| Spearman p < 0.05 | 65 |
| Theil–Sen 95 % CI entirely below 0 | 56 (none entirely above) |

(The metrics table rounds p to three decimals; recomputed from the cube the BH count is 47, from the rounded
table 48. Quote 47.)

**What it says.** About half the ranges show an individually significant relationship, and 47 of them
survive a correction for having looked 160 times. Under the null, ~8 would.

**Significance and caveat.** FDR handles *multiplicity*. It does not handle *dependence*: the ranges share
the same eleven years and neighbouring ranges share weather (section 8), so the 160 tests are not 160
independent looks. Treat "74 of 160" and "47 after FDR" as descriptive counts, and let section 4 carry the
inferential weight.

**Why report it.** One sentence in Methods or the supplement. It answers the reviewer who asks "how many of
these would be significant by chance?" without being the main evidence.

---

## 3. The robust pair: Theil–Sen slope with CI, Spearman ρ

**How.** Per range, the Theil–Sen slope (median of the 55 pairwise slopes for n = 11) with its
distribution-free 95 % CI, and Spearman's rank correlation with its p-value. Computed in `per_range_fits`
above; stored in the metrics table as `theil_sen_slope`, `theil_sen_low`, `theil_sen_high`, `spearman_rho`,
`spearman_pval`; drawn by the robust twin of the world map (`plot_temperature_sensitivity_map(robust=True)`).

**Why.** With eleven points, one extreme year (WY2015 in western North America; a −36-day Sierra Nevada
year) can pull an OLS line and inflate r. Rank-based and median-based statistics are insensitive to that,
so agreement between the two families is the cheapest available evidence that the pattern is not an artefact
of a few points.

**What it says.**

| comparison over the 160 ranges | value |
| --- | --- |
| Theil–Sen slope: 5th / median / 95th percentile | −12.4 / **−5.0** / +3.3 days/°C |
| share with negative Theil–Sen slope | 0.83 |
| median \|OLS − Theil–Sen\| slope difference | 0.9 days/°C (90th percentile 3.6) |
| ranges where the two slopes differ in sign | 11 — all with \|r\| < 0.55: the three Central Andes cordilleras, Hawaii, Kaikoura, Daba, Qin Ling, Alashan, Stanovoy, Northern Andes |
| median Spearman ρ | −0.57 (share negative 0.83) |
| ranges with \|r\| ≥ 0.6 but \|ρ\| < 0.4 (one-point-driven Pearson) | 2 |
| ranges with \|ρ\| ≥ 0.6 but \|r\| < 0.4 | 5 |
| median half-width of the Theil–Sen 95 % CI | 7.5 days/°C |

**Significance.** Wherever the relationship is worth talking about, the estimator does not matter. The
Theil–Sen CI half-width of 7.5 days/°C is the honest statement of how precisely *one* range's sensitivity is
known from eleven years; it is why per-range values should be presented with r and n rather than as
standalone measurements, and why the CI itself should not be used to colour or threshold a map (it is built
from 55 pairwise slopes and is very coarse).

**Why report it.** One sentence in Methods ("Theil–Sen slopes agree with OLS within 1 day/°C for half the
ranges and change sign only in 11 weakly correlated ranges") plus the robust twin map and the columns in the
supplementary table. It pre-empts the "outlier year" critique.

---

## 4. Field significance: the shared-year permutation test (the decisive test)

**How.** Shuffle the water-year labels of the *temperature* series, using **one permutation applied
identically to every range**, recompute all 160 fits, and record four summary statistics: the share of
ranges with a negative slope, the median slope, the count with p < 0.05, and the pooled slope (section 5).
Repeat 2000 times. The p-value is the share of shuffles at least as extreme as the observed value.

**Why this design.** The ranges are not independent: they share the same eleven years, and neighbouring
ranges see the same weather (section 8). Shuffling years the same way for every range keeps that entire
dependence structure intact: whatever regional coherence exists in the real data exists in every shuffle.
Only the pairing between temperature and onset is broken. The test therefore answers exactly the question a
sceptical reader asks: *could this many negative slopes arise by chance from eleven shared years?*

```python
rng = np.random.default_rng(0)
temp_by_range = pairs.pivot(index='range', columns='water_year', values='x')
pooled = sps.linregress(pairs['x'], pairs['y'])
observed = np.array([(fits.slope < 0).mean(), fits.slope.median(), (fits.p < 0.05).sum(), pooled.slope])
null = []
for _ in range(2000):
    perm = dict(zip(years, rng.permutation(years)))                          # ONE year relabelling for all ranges
    shuffled = pairs.copy()
    shuffled['x'] = temp_by_range.to_numpy()[temp_by_range.index.get_indexer(shuffled['range']),
                                            temp_by_range.columns.get_indexer(shuffled['water_year'].map(perm))]
    shuffled = shuffled.dropna(subset=['x'])
    f = per_range_fits(shuffled)
    null.append([(f.slope < 0).mean(), f.slope.median(), (f.p < 0.05).sum(),
                 sps.linregress(shuffled['x'], shuffled['y']).slope])
null = np.array(null)
p_values = [(null[:, 0] >= observed[0]).mean(), (null[:, 1] <= observed[1]).mean(),
            (null[:, 2] >= observed[2]).mean(), (null[:, 3] <= observed[3]).mean()]
```

**What it says (v10, 2000 permutations).**

| statistic | observed | null 95th (or 5th) percentile | most extreme null | permutation p |
| --- | --- | --- | --- | --- |
| share of ranges with slope < 0 | 0.86 | 0.65 | 0.79 | < 0.0005 |
| median per-range slope | −4.9 days/°C | −1.3 | −2.6 | < 0.0005 |
| ranges with p < 0.05 | 74 | 15 | 30 | < 0.0005 |
| pooled slope | −4.6 days/°C | −1.3 | −2.9 | < 0.0005 |

No shuffle in 2000 came close on any statistic ("< 0.0005" is the resolution of 2000 draws, 1/2000).

Note what the null *does* produce: up to 30 "significant" ranges and up to 79 % negative slopes purely from
eleven shared years and regional coherence. That is why counting significant ranges (section 2) is not
evidence on its own.

**Significance.** This is the statement that "warmer springs bring earlier onset across the world's
mountains" rests on. It is valid with eleven years because the inference is across 160 ranges under a null
that preserves their dependence, not across eleven years in one range.

**Limits.** It tests the global null of *no relationship anywhere*. It does not say which ranges are real
(section 2 and the per-range table do that, descriptively), and with only 11! ≈ 4×10⁷ distinct year
permutations available, 2000 draws is a sample of them (adequate; the observed values are far outside the
null).

**Why report it.** Main text, one sentence with the p-value; Methods, one paragraph describing the design;
supplement, a histogram of the null for each statistic with the observed value marked.

---

## 5. One global effect size: the pooled slope with a year-block bootstrap CI

**How.** OLS through all 1580 range-years at once. Because both variables are already anomalies about each
range's own median, this is a within-range (fixed-effects) regression: it asks, within a range, how much
earlier onset is in a warmer-than-usual spring, averaged over all ranges. The CI comes from a block bootstrap
by water year: resample the eleven years with replacement, take *all* ranges' values for the drawn years,
refit, repeat 2000 times, and read the 2.5th and 97.5th percentiles.

**Why.** Any write-up needs one number with an honest interval. The naive OLS standard error treats the 1580
points as independent; they are not (160 ranges × 11 shared years). Resampling *years* keeps every range's
values for a year together, so the interval reflects the eleven independent things we actually have.

```python
t = sps.t.ppf(0.975, len(pairs) - 2)
naive_ci = (pooled.slope - t * pooled.stderr, pooled.slope + t * pooled.stderr)

def year_block_bootstrap(stat, n_boot=2000):
    out = []
    for _ in range(n_boot):
        draw = rng.choice(years, len(years), replace=True)                    # resample WATER YEARS
        out.append(stat(pd.concat([pairs[pairs.water_year == y] for y in draw])))
    return np.percentile(out, [2.5, 97.5])

pooled_ci = year_block_bootstrap(lambda b: sps.linregress(b['x'], b['y']).slope)
median_slope_ci = year_block_bootstrap(
    lambda b: per_range_fits(b.drop_duplicates(['range', 'water_year'])).slope.median(), n_boot=1000)
```

**What it says (v10).**

| | value |
| --- | --- |
| pooled slope | **−4.6 days/°C** |
| pooled Pearson r (1580 range-years) | −0.41 (r² = 0.17) |
| naive 95 % CI (independent points assumed) | −5.1 to −4.1 (±0.5) |
| **year-block bootstrap 95 % CI** | **−5.5 to −3.7** (±0.9) |
| median per-range slope, year-block bootstrap 95 % CI | −4.9; −5.7 to −4.4 |

The naive interval is 1.8× too narrow. The pooled r of −0.41 says spring temperature explains about a sixth
of the range-year variance globally: a moderate control at the global scale, much stronger in the
mid-latitudes (section 7), weaker at high latitudes and in the tropics.

**Significance.** Globally, a spring 1 °C warmer than a range's usual brings onset 4 to 5 days earlier, with
an interval that honestly reflects eleven years. It is also the number that lets the study be compared with
station-based sensitivities in the literature, provided the comparison notes that this is *interannual*
sensitivity to spring temperature, not a trend and not a warming response (section 12).

**Why report it.** Main text, with the bootstrap interval. Never quote the naive interval.

---

## 6. Leave-one-year-out

**How.** Drop one water year, recompute the pooled slope, the median per-range slope, the share of negative
slopes and the count of significant ranges; repeat for each of the eleven years.

**Why.** WY2015 was an extreme early-melt year in western North America (the cross-range mean onset anomaly
that year is −8 days; no other year exceeds −3), and a reviewer will ask whether it makes the result.

```python
for y in years:
    sub = pairs[pairs.water_year != y]
    f = per_range_fits(sub); o = sps.linregress(sub['x'], sub['y'])
    print(y, round(o.slope, 2), round(f.slope.median(), 2), round((f.slope < 0).mean(), 2), (f.p < 0.05).sum())
```

| dropped year | pooled slope | median per-range slope | share negative | ranges p < 0.05 |
| --- | --- | --- | --- | --- |
| WY2015 | −4.57 | −4.88 | 0.88 | 71 |
| WY2016 | −4.34 | −4.95 | 0.86 | 67 |
| WY2017 | −4.83 | −4.94 | 0.85 | 71 |
| WY2018 | −4.74 | −5.04 | 0.86 | 69 |
| WY2019 | −4.58 | −4.71 | 0.83 | 66 |
| WY2020 | −4.52 | −4.84 | 0.82 | 66 |
| WY2021 | −4.80 | −4.99 | 0.82 | 68 |
| WY2022 | −4.45 | −4.91 | 0.84 | 73 |
| WY2023 | −4.39 | −4.81 | 0.84 | 70 |
| WY2024 | −4.61 | −5.15 | 0.87 | 66 |
| WY2025 | −4.63 | −4.90 | 0.85 | 69 |

**What it says.** Nothing depends on any one year. Dropping WY2015 leaves the pooled slope at −4.57 days/°C
and the share of negative slopes at 0.88.

**Why report it.** One sentence in Methods ("the pooled slope stays between −4.3 and −4.8 days/°C when any
single water year is withheld"), the table in the supplement.

---

## 7. Structure with latitude

**How.** Group ranges by absolute centroid latitude; compare per-range slopes between the mid-latitude and
polar groups with a Mann–Whitney U test (rank-based, no normality assumption); Spearman correlation of
|latitude| with slope over all ranges.

**Why.** The physical story is that mid-latitude melt is paced by spring air temperature while high-latitude
melt near the solstice is paced by insolation. The latitude structure of the sensitivity is that claim's
statistical form.

```python
fits['abslat'] = fits['lat'].abs()
bands = {(0, 35): 'tropics and subtropics', (35, 50): 'mid-latitudes', (50, 60): 'boreal', (60, 90): 'polar'}
for (lo, hi), name in bands.items():
    b = fits[(fits.abslat >= lo) & (fits.abslat < hi)]
    print(name, len(b), b.slope.median(), b.r.median(), (b.p < 0.05).mean())
mid, pol = fits[(fits.abslat >= 35) & (fits.abslat < 50)], fits[fits.abslat >= 60]
mw = sps.mannwhitneyu(mid.slope, pol.slope, alternative='less')
sp = sps.spearmanr(fits.abslat, fits.slope)
```

| \|latitude\| band | ranges | median slope [days/°C] | median r | share p < 0.05 |
| --- | --- | --- | --- | --- |
| 0–35° | 30 | −4.8 | −0.48 | 0.33 |
| 35–50° | 76 | **−6.8** | **−0.67** | **0.62** |
| 50–60° | 26 | −4.6 | −0.63 | 0.38 |
| 60–90° | 28 | −2.4 | −0.36 | 0.25 |

Mann–Whitney, mid-latitude slopes more negative than polar: p = 2 × 10⁻⁶. Spearman of |latitude| against
slope over all 160 ranges: ρ = 0.24, p = 0.002 (positive: slopes get less negative poleward).

**What it says.** The sensitivity is strongest and most consistent between 35° and 50°, and about a third as
strong poleward of 60°. The tropical group is intermediate in slope but weakly correlated: high variability,
few of the fits significant.

**Why report it.** Main text, as the physical result, with the Mann–Whitney p-value. It is a between-group
comparison of 76 vs 28 ranges, so the small-n-per-range problem does not touch it.

---

## 8. How dependent are the ranges? (diagnostics that justify sections 4 and 5)

**How.** Pairwise correlations of the spring-temperature anomaly series between ranges; the share of onset
and temperature variance that is common to all ranges in a year; the cross-range mean onset anomaly per
year; and, for contrast, the sign test that a naive analysis would use.

```python
corr_t = temp_by_range.T.corr(min_periods=6).values
iu = np.triu_indices_from(corr_t, 1)
mean_pairwise, share_strong = np.nanmean(corr_t[iu]), (np.abs(corr_t[iu]) > 0.5).mean()
onset_by_range = pairs.pivot(index='range', columns='water_year', values='y')
common_share = 1 - np.nanvar(onset_by_range.sub(onset_by_range.mean(axis=0), axis=1).values) / np.nanvar(onset_by_range.values)
yearly_mean = onset_by_range.mean(axis=0)
k = int((fits.slope < 0).sum())
sign_test_p = sps.binomtest(k, len(fits), 0.5, alternative='greater').pvalue   # assumes independent ranges: too small
```

| diagnostic | value |
| --- | --- |
| mean pairwise inter-range correlation of the spring-T anomaly | 0.04 |
| share of range pairs with \|correlation\| > 0.5 | 0.20 |
| share of onset-anomaly variance in the all-range yearly mean | 0.03 |
| share of spring-T variance in the all-range yearly mean | 0.05 |
| cross-range mean onset anomaly, WY2015 … WY2025 [days] | −8.1, −1.5, −0.9, −0.8, −0.2, −3.0, −2.9, −1.9, −0.3, −1.8, −0.5 |
| binomial sign test on 137/160 negative slopes (wrongly assuming independence) | p = 3 × 10⁻²¹ |

**What it says.** Dependence is *regional*, not global: one range pair in five is strongly correlated
(western North America, High Mountain Asia, the Alps–Carpathians …), but no single global "year signal"
dominates (3–5 % of variance). So the 160 ranges are neither independent nor one big replicate; they are
perhaps a few dozen effectively independent clusters. That is precisely the situation the shared-year
permutation (section 4) and the year-block bootstrap (section 5) are built for, and precisely why the sign
test's p = 10⁻²¹ and the naive CI must not be quoted.

**Why report it.** Methods, briefly, as the justification of the test design; the pairwise-correlation
statement is enough.

---

## 9. The worked example (Sierra Nevada)

| statistic | value |
| --- | --- |
| n | 11 |
| OLS slope, r, p | −12.6 days/°C, −0.84, 0.001 |
| Theil–Sen slope [95 % CI] | −11.8 [−16.8, −4.6] |
| Spearman ρ | −0.73 |

One of the few ranges whose relationship is individually significant beyond doubt, with OLS and Theil–Sen in
agreement and a CI that, for once, is informative. It is the right exemplar because its outlier year (WY2015,
−36 days) is *visible* in the scatter and the robust slope shows it is not doing the work.

---

## 10. OLS or robust: what to show and what to report

**Show OLS on the main map; report both; put the robust twin in the supplement.**

- OLS slope and Pearson r are what readers expect, are the statistics whose p-values feed sections 2 and 4,
  and match the worked example. The robust twin is visually almost identical for every range that carries the
  story (section 3), so showing it as the main figure would buy nothing while making the caption harder.
- The robust statistics do their work as a *check*: one Methods sentence, the supplementary figure, and two
  columns in the supplementary table. Their existence is the answer to "isn't this driven by 2015?".
- Do not report R² per range (with n = 11 it is dominated by sampling noise), do not report the Theil–Sen CI
  as a significance criterion (56 ranges exclude zero, but the interval is built from 55 pairwise slopes), and
  do not colour any map by p-value.
- If a single robust number is wanted in the main text, the median Theil–Sen slope across ranges (−5.0
  days/°C) sits next to the OLS median (−4.9) and the pooled slope (−4.6): three routes, one answer.

---

## 11. Suggested reporting

**Main text (three sentences).** *Across 160 mountain ranges with at least three years of paired data, the
annual runoff onset anomaly decreased with the spring temperature anomaly in 86 % of ranges (median OLS
slope −4.9 days/°C, IQR −7.4 to −2.0; median r = −0.60). A permutation test that reassigns water-year labels
identically across all ranges, preserving their spatial dependence, never reproduced the observed fraction of
negative slopes, median slope or number of significant ranges in 2000 draws (p < 0.0005); a pooled
within-range regression gives −4.6 days/°C (year-block bootstrap 95 % CI −5.5 to −3.7), insensitive to
withholding any single year. Sensitivity peaks between 35° and 50° latitude (median −6.8 days/°C, 62 % of
ranges p < 0.05) and weakens poleward of 60° (−2.4 days/°C; Mann–Whitney p < 10⁻⁵).*

**Methods (one paragraph).** The pairs and their rules (section 0); per-range OLS with r, p, n, and the
Theil–Sen / Spearman check; BH-FDR count; the shared-year permutation design with its four statistics; the
year-block bootstrap; leave-one-year-out; the latitude tests. State that "sensitivity" means interannual
sensitivity to spring temperature.

**Supplement.** The per-range table (slope, r, n, p, q, Theil–Sen slope and CI, ρ); the robust twin map;
the leave-one-year-out table; histograms of the four permutation nulls with the observed values marked; the
histogram of per-range slopes and r.

---

## 12. Limitations and anticipated critiques

1. **Eleven years.** Every per-range statistic is imprecise (slope ±7 days/°C, r ±0.2). *Response:* the
   inference is across 160 ranges under a dependence-preserving null; per-range values are presented as
   descriptive with n and r; the worked example is one of the few individually solid ranges.
2. **The anomaly baseline is the same eleven years.** Each anomaly is relative to the median of the same
   record, which removes a degree of freedom and makes the anomalies slightly negatively dependent within a
   range. *Response:* this affects all ranges equally, is inside the permutation null, and the median (not
   mean) baseline limits the effect of extreme years.
3. **Spatial dependence.** Neighbouring ranges share weather, so "74 of 160" is inflated. *Response:*
   section 8 quantifies it, section 4 is designed for it, section 5's interval respects it.
4. **Correlation is not mechanism.** Spring temperature covaries with cloud cover, humidity, longwave
   radiation and precipitation; a warm spring is also a sunny or a dry one in many climates. *Response:* the
   statistic is a *sensitivity*, an empirical slope, not an attribution; the regional pattern (strong where
   melt is thermally paced, weak where it is insolation-paced) is consistent with mechanism but does not prove
   it. The monthly regressions on the other ERA5-Land variables exist for a fuller picture.
5. **Not a warming sensitivity and not a trend.** Interannual sensitivity to spring temperature is not the
   same quantity as the response of melt timing to long-term warming (which also shifts accumulation and the
   annual cycle), and eleven years support no trend statement at all. *Response:* say "interannual" and never
   "trend".
6. **ERA5-Land temperature over complex terrain.** A 9 km reanalysis smooths mountain temperature; the
   anomaly (not the absolute) is used, which is far less biased, but lapse-rate and inversion effects are
   not resolved. *Response:* state it; the per-range averaging over the range's snow-covered pixels reduces
   it further.
7. **Selection of ranges and years.** Range-years with too few valid pixels are dropped; ranges with < 3
   pairs are excluded; the 30 ranges with insets are a curated, legibility-driven subset. *Response:* the
   rules are stated and identical to the metrics table; the statistics use all 160 ranges, not the 30.
8. **Choice of spring months.** MAM / SON is fixed; for early-melting ranges the relevant window may be
   earlier. *Response:* the month-by-month regressions show which months carry the signal per range; the
   spring mean is a compromise chosen a priori, not tuned.
9. **Positive sensitivities.** About 14 % of slopes are positive, one significantly (Verkhoyansk). *Response:*
   nearly all are in weakly correlated ranges where the sign is unresolved; the Verkhoyansk case belongs in
   the discussion as a high-latitude range where warm springs may mean more snowfall or later snow-off
   rather than as a counterexample.
10. **Permutation resolution.** 2000 draws cannot resolve p below 1/2000. *Response:* the observed values are
    far outside the null; increasing the draws would only shrink the reported bound.

---

## 13. Housekeeping

- Sections 4–8 should be ported into `temperature_sensitivity.ipynb` (the code above runs as-is after the
  notebook's setup cells) and their tables written to `results/<version>/`, so that the published numbers are
  regenerated with every dataset version like the metrics table is.
- The metrics table rounds to three decimals; recompute q-values from the cube, not from the table.
- The colorbar of the world maps shows the slope colour but not the transparency that encodes r; the
  caption must say so.
- Numbers in this document are for v10. Rerun the blocks when the dataset version changes and update the
  tables; the v10 values are also recorded with dates in this file's history.
