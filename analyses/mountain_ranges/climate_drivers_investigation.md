# Beyond spring temperature: what else paces the interannual runoff-onset anomaly, and how to find out properly

Status: investigation note written 2026-09-08 against dataset version **v10** (water years 2015–2025). Every number
comes from the scripts in `scratch/climate_drivers/` (gitignored; `README.md` there says how to rerun them), run on
the v10 mountain-range cube (`data/aggregation/v10/all_mountain_ranges_fcf_lte_50.nc`) under the same "analyses'
view" rules as `temperature_sensitivity.ipynb` (bins ≤ 100 pixels masked, a bin-year needs > 30 % of the bin's median
pixels, the tropical-Andes rule, a range-year mean anomaly needs ≥ 10 % of the range's median pixels). Nothing here is
in a notebook yet. Companion to `spring_temperature_sensitivity_stats.md`, which establishes *that* warmer springs bring
earlier onset; this note asks *what else does, where, and when*, and ends with a plan for doing that methodically.

The prompt: spring 2 m temperature explains most of the year-to-year onset variation in ranges like the Sierra Nevada
and almost none of it elsewhere. Does that mean other variables pace the onset elsewhere? Which, in which months, how
much variance, are they independent of temperature, and what is missing from the eight ERA5-Land variables?

---

## Summary

1. **The eleven ERA5-Land predictors are two signals, not eleven.** Pooled across all range-years, a principal-component
   analysis of the eight spring-mean anomalies puts 70 % of their variance on two axes: a *warm–humid* axis (2 m
   temperature, dewpoint, downwelling longwave together, opposite the latent heat flux) and a *snowy–cloudy vs
   sunny* axis (snow water equivalent, snowfall, sensible heat flux, longwave together, opposite downwelling
   shortwave). Dewpoint and longwave are near-duplicates of temperature (pooled r = 0.65 and 0.57 with T, 0.87 with
   each other). The three "energy sum" variables add no information and one of them is mis-specified (section 2).
2. **Spring temperature and the snow amount are the two real drivers; downwelling radiation anomalies are not.**
   Median r² across ranges at the spring mean: temperature 0.40, snow water equivalent 0.35, latent heat flux 0.32,
   snowfall 0.24, sensible heat 0.17, dewpoint 0.17, longwave 0.12, **shortwave 0.06** (the no-relationship level for
   eleven years is about 0.10). Monthly downwelling-shortwave anomalies explain nothing anywhere, in any month or
   latitude band. The interannual timing is not paced by how sunny the spring was; insolation sets the *mean*
   timing (the latitude/aspect results), not the anomaly.
3. **The turbulent fluxes and the spring-mean SWE are consequences of the onset, not causes.** ERA5 fluxes are positive
   downward. An early melt-out exposes ground that heats and dries the air, so the spring sensible and latent fluxes go
   strongly negative (upward) and the spring SWE is low; a late onset leaves the surface snow-covered and both fluxes
   positive. That is exactly the sign pattern observed (r = +0.29, +0.55, +0.58). Read them as diagnostics of the snow
   state during the window, not as controls, and do not put them on the predictor side of a regression.
4. **Variable importance shifts with latitude far more than with month.** Pooled within-range R² by |latitude| band:

   | band | ranges | spring T | March SWE | spring snowfall | spring SWdown | winter T | T + SWE + snowfall |
   | --- | --- | --- | --- | --- | --- | --- | --- |
   | 0–35° | 30 | 0.14 | 0.07 | **0.18** | 0.00 | 0.01 | 0.24 |
   | 35–50° | 76 | **0.31** | 0.23 | 0.16 | 0.01 | 0.07 | 0.39 |
   | 50–60° | 26 | **0.31** | 0.17 | 0.22 | 0.03 | 0.19 | 0.47 |
   | 60–90° | 28 | 0.07 | **0.34** | 0.17 | 0.01 | 0.00 | 0.46 |

   Poleward of 60° the snow amount at the start of spring (standardised coefficient +0.58) outweighs temperature
   (−0.26) by a wide margin: deeper snow, later runoff onset, whatever the spring temperature. In the tropics and
   subtropics spring snowfall is the only driver with any signal. Between 35° and 60° temperature leads and the snow
   amount adds about a quarter more explained variance.
5. **Within the months, the signal sits in February–April (NH), and the three-month spring mean beats any single
   month.** Temperature's median r² by month is 0.06 (Dec), 0.08, 0.23 (Feb), 0.22 (Mar), 0.19 (Apr), 0.08 (May) and
   0.40 for the spring mean. Re-centring the window on each range's own onset month does not help on average (0.35
   for the onset month alone); it helps only the 27 ranges whose median onset falls before March (0.19 → 0.26).
6. **How much variance can be explained in a given range: less than the heatmaps suggest.** The best of the 99
   (variable × month) cells reaches a median r² of 0.70 across ranges, but the same search on onset series with the
   years shuffled reaches **0.52**. Per-range "which variable wins" is mostly selection noise with eleven years. Under a
   pre-specified five-driver set the observed median best r² is 0.52 against a shuffled 0.21, which is a real signal,
   but the honest per-range figure is the leave-one-year-out predictive R²: median 0.02 for spring temperature alone
   (0.38 in the 47 ranges with r² ≥ 0.5). Per-range multivariate fits are not supportable at n = 11; pooled and
   hierarchical models are (section 6).
7. **Temperature and snow amount are confounded in most mid-latitude ranges, including the Sierra Nevada.** Within
   range, spring T and March SWE correlate at a median −0.51 (0–35°), −0.44 (35–50°), −0.36 (50–60°) and −0.06
   (60–90°). In the Sierra Nevada r(T, SWE) = −0.77: its 0.71 r² for temperature is 0.63 for snow amount, and the two
   cannot be separated with these eleven years (WY2015 was both the warmest and the least snowy). The Cascades
   (r = −0.34) separate them and both matter. Of the pooled R² of 0.35 for T + SWE + snowfall, 0.18 is shared.
8. **For about a third of the weakly temperature-sensitive ranges, nothing in ERA5-Land explains the onset anomaly.**
   The polar and Kamchatka group is explained by SWE (section 5). The dry interior-Asian plateaus, the tropical Andes,
   the New Zealand and Tasmanian ranges, and the low, marginal-snow European ranges show r² ≤ 0.2 for every a priori
   driver. That points at a noise floor (Sentinel-1 revisit, few pixels, the anomaly baseline), at 9 km reanalysis not
   representing those ranges, or at processes the eight variables do not carry (section 7), rather than at another
   ERA5 variable.

---

## 1. The data and the analyses' view (unchanged from the notebooks)

```python
# scratch/climate_drivers/build_pairs.py, in short
cube = xr.open_dataset(paths.aggregation_dir('mountain_ranges', 'v10') / 'all_mountain_ranges_fcf_lte_50.nc')
ds = aggregate.collapse(cube); ds = aggregate.threshold(ds, 101)                 # CHILI folded, bins <= 100 px masked
# ... bin-year > 30 % rule, tropical-Andes rule (as in temperature_sensitivity.ipynb cell 4) ...
onset_anomaly = aggregate.weighted_mean(ds, 'runoff_onset_anomaly', ['elevation', 'aspect']).where(valid_fraction >= 0.1)
climate = ds[era5.VARIABLES]                                                      # (range, water_year, month), 8 anomalies
```

161 ranges, 1589 range-years with an onset anomaly (the stats note's 160/1580 after its ≥ 3-pairs rule). The cube's
month axis is hemisphere-aware: `winter_month_1..3` = Dec–Feb (NH) / Jun–Aug (SH), `spring_month_1..3` = Mar–May /
Sep–Nov, `summer_month_1..3` = Jun–Aug / Dec–Feb. Every ERA5-Land value is an anomaly (the cell minus its per-pixel
median over the eleven water years), zonally averaged over the range with weights = polygon coverage × seasonal-snow
fraction × onset validity (`era5.zonal_anomalies`).

---

## 2. What the variables are, physically, and what role each can play

ERA5-Land is ECMWF's land-surface reanalysis at 0.1° (~9 km): an offline run of the HTESSEL land model forced by ERA5
atmosphere. **Temperature, dewpoint, radiation and snowfall are atmospheric forcing** (from ERA5, which assimilates
observations). **Snow water equivalent, the turbulent fluxes, and everything about the snowpack are model output** of a
single-layer snow scheme with no data assimilation over land. The ECMWF sign convention is *positive downward* for all
surface fluxes; the `_sum` variables are monthly sums of daily totals in J m⁻² (radiation, fluxes) or m of water
(snowfall); the others are monthly means.

| variable | what it is | what a positive spring anomaly means for a snowpack | a priori role | median r (spring mean) |
| --- | --- | --- | --- | --- |
| `temperature_2m` | 2 m air temperature, monthly mean [K] | warmer air: more sensible heat to the snow, higher LWdown, earlier ripening; rain instead of snow | **driver** | −0.62 |
| `dewpoint_temperature_2m` | 2 m dewpoint [K]: absolute humidity | moister air: condensation on the snow (large latent release), higher LWdown; nearly the same information as T (r = 0.65) plus cloud/humidity | driver, redundant with T | −0.35 |
| `surface_thermal_radiation_downwards_sum` | downwelling longwave [J m⁻²]: emitted by clouds, water vapour, warm air | cloudier / warmer / moister sky; the largest energy input to snow under cloud | driver; tracks Td (0.87) and T (0.57) | −0.27 |
| `surface_solar_radiation_downwards_sum` | downwelling shortwave [J m⁻²] | sunnier spring; but snow reflects 60–90 % of it, so the *absorbed* fraction depends on albedo, which is not in the store | driver in principle; **no signal in practice** | +0.02 |
| `snowfall_sum` | snowfall, monthly total [m w.e.] | more snow falling in the window: fresh high-albedo surface, added cold content, a wetted surface re-frozen; also a colder, cloudier month (r with T −0.38) | driver (late-season) | +0.48 |
| `snow_depth_water_equivalent` | modelled SWE, monthly mean [m w.e.] | more snow on the ground; **in March (NH): the accumulation-season outcome, an antecedent state**; in April–May: partly the *result* of melt timing | **driver as `spring_month_1`**; state/consequence as the spring mean | +0.58 (spring mean) |
| `surface_sensible_heat_flux_sum` | turbulent heat, air → surface positive [J m⁻²] | more heat delivered to the surface; but the spring-mean anomaly is dominated by whether the ground was snow-covered (SH > 0 over melting snow, strongly < 0 over bare warm ground) | **consequence** of onset timing | +0.29 |
| `surface_latent_heat_flux_sum` | turbulent moisture flux, condensation positive [J m⁻²] | condensation onto the snow (energy in) vs sublimation/evapotranspiration (energy out); after melt-out, evapotranspiration makes it strongly negative | **consequence** | +0.55 |
| `total_radiation_sum` (derived) | SWdown + LWdown | mostly the SWdown anomaly by magnitude; SW and LW anomalies anticorrelate (−0.66, clouds), so the sum cancels the cloud signal | drop | +0.07 |
| `total_turbulent_flux_sum` (derived) | SH + LH | a snow-cover indicator (see SH, LH) | consequence | +0.46 |
| `total_surface_energy_sum` (derived) | SWdown + LWdown + SH + LH | **not a surface energy balance**: it omits reflected shortwave (the largest term over snow) and emitted longwave (~−300 W m⁻² at 0 °C), so it is dominated by the term the snow mostly reflects | drop; use ERA5-Land's `surface_net_solar_radiation_sum` and `surface_net_thermal_radiation_sum` instead | +0.27 |

Two consequences for the existing heatmaps. First, the eleven figures show two or three independent things: the
temperature/dewpoint/longwave figures are one story, the SWE/snowfall figures another, and the four flux and sum
figures mostly show the onset anomaly correlated with itself through the snow-cover state. Second, the heatmap of
`total_surface_energy_sum` should not be interpreted as "the energy available for melt".

---

## 3. Variable × month screen: where the signal is (and what the null looks like)

**How.** For every range, variable and month (the nine months and the three seasonal means), OLS of the range mean
onset anomaly on the anomaly, n ≥ 6 (as the notebook since 2026-09-05). Summaries across the 151 ranges with fits.

**Median r² across ranges** (rows: variable; columns: NH-equivalent month; the last three are the seasonal means):

| | Dec | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | winter | **spring** | summer |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T2m | .06 | .08 | **.23** | **.22** | .19 | .08 | .12 | .05 | .08 | .14 | **.40** | .07 |
| Td2m | .08 | .08 | .13 | .11 | .12 | .08 | .09 | .06 | .09 | .11 | .17 | .09 |
| SWE | .10 | .12 | **.23** | **.29** | **.30** | **.28** | .18 | .15 | .06 | .16 | **.35** | .19 |
| snowfall | .07 | .08 | .14 | .15 | .11 | .09 | .06 | .09 | .05 | .20 | .24 | .07 |
| SWdown | .08 | .05 | .07 | .08 | .07 | .09 | .07 | .06 | .06 | .08 | .06 | .05 |
| LWdown | .09 | .05 | .08 | .09 | .08 | .08 | .08 | .05 | .10 | .08 | .12 | .07 |
| SH | .07 | .06 | .09 | .15 | .14 | .13 | .12 | .08 | .04 | .08 | .17 | .09 |
| LH | .07 | .06 | .11 | .18 | .21 | .14 | .12 | .08 | .06 | .09 | .32 | .09 |

The reference level: a single fixed predictor with no relationship gives a median r² of about 1/(n−1) ≈ 0.10 at
n = 11. Everything at 0.05–0.10 is indistinguishable from noise: shortwave and longwave in every month, dewpoint in most
months, and every variable in December–January and July–August. The share of ranges with p < 0.05 tells the same story
(chance = 0.05): spring-mean T 0.48, SWE 0.46, LH 0.42, snowfall 0.26, Td 0.26, SH 0.28, LWdown 0.13, SWdown 0.07.

**The sign is consistent where the signal is.** Spring-mean T negative in 88 % of ranges, SWE positive in 93 %,
snowfall positive in 88 %, LH positive in 81 %; SWdown is negative in 48 %, i.e. a coin flip.

**The selection null.** For each range take the best of the 8 × 9 single-month cells:

| | observed | onset years shuffled within range (20 repeats) |
| --- | --- | --- |
| median over ranges of the best-of-99 r² | 0.70 | 0.52 ± 0.02 |
| median over ranges of the best-of-5 a priori drivers (fixed windows) | 0.52 | 0.21 ± 0.02 |
| median spring-T r² (one fixed cell) | 0.40 | ~0.10 |

So a heatmap cell with r² ≈ 0.6 in a single month for a single range is what noise produces once you have looked at 99
cells. The month structure is real *in aggregate* (Feb–Apr for T; Feb–May for SWE), but per range the "best month"
should not be read off a heatmap. Where the best cell sits, in aggregate: Apr 35 ranges, May 27, Mar 26, Jun 21, Feb 16;
July–August and December–January together 26.

```python
# the selection null, in short (scratch/climate_drivers/a_univariate.py)
shuffled['onset_anomaly'] = shuffled.groupby(['mountain_range', 'month'])['onset_anomaly'].transform(
    lambda s: s.sample(frac=1, random_state=rep).values)          # break the pairing, keep everything else
best_r2 = max(np.corrcoef(g[v], g['onset_anomaly'])[0, 1] ** 2 for v in VARS for month ...)   # per range
```

---

## 4. The month question, and the onset-relative window

**How.** Each range's median onset (pixel-weighted over its bins) converted to a calendar month; predictors taken at
that month (m0), the month before, two before, and their means; compared with the fixed spring mean.

| median onset month (NH-equivalent) | Jan | Feb | Mar | Apr | May | Jun |
| --- | --- | --- | --- | --- | --- | --- |
| ranges | 5 | 25 | 61 | 47 | 21 | 2 |

129 of 161 ranges melt inside the fixed spring window; 30 melt before it (Jan–Feb: the Alashan Plateau, Qin Ling, the
Adirondacks, the Appalachians, Hokkaido, the Southern Alps…).

| median r² across ranges | fixed spring mean | m0 − 2 | m0 − 1 | m0 | mean(m0−1, m0) | mean(m0−2 … m0) |
| --- | --- | --- | --- | --- | --- | --- |
| T2m | **0.40** | 0.06 | 0.17 | 0.35 | 0.31 | 0.27 |
| SWE | **0.35** | 0.14 | 0.14 | 0.35 | 0.26 | 0.18 |
| snowfall | 0.24 | 0.08 | 0.14 | 0.17 | 0.22 | 0.24 |
| LWdown | 0.12 | 0.05 | 0.07 | 0.14 | 0.14 | 0.10 |

| T2m, by where the range's onset falls | spring mean | m0 − 1 | m0 | mean(m0−1, m0) | ranges |
| --- | --- | --- | --- | --- | --- |
| before the spring window | 0.19 | 0.09 | **0.26** | 0.15 | 27 |
| inside it | **0.44** | 0.19 | 0.40 | 0.36 | 122 |

**What it says.** The signal is concentrated in the onset month itself and the fixed three-month mean is a good
estimator of it: averaging three months removes month-to-month noise and the median range melts in March–April
anyway. The a priori spring window is not what makes 50 ranges look insensitive. For the 27 early-melting ranges an
onset-centred window recovers some signal (0.19 → 0.26) and should be used in a follow-up, but it does not turn them into
Sierra Nevadas. Note also that the month *before* onset carries little (0.17): the anomaly is set close to the event,
consistent with a threshold process (the first sustained ripening) rather than a cumulative one.

**Where winter matters.** In the 50–60° band the winter-mean temperature and dewpoint have r² 0.35 and 0.37 with
negative r (−0.59, −0.61): warm, moist boreal winters go with early onset, most plausibly through mid-winter melt and
rain-on-snow events and through a shallower, warmer pack. Poleward of 60° the winter SWE and snowfall dominate instead
(r² 0.38, 0.30; positive). In the 0–35° band the winter T anomaly has a *positive* median r (+0.23): warm winters, later
onset, which needs an explanation (a monsoon/precipitation covariate, or which pixels survive a warm winter).

---

## 5. Co-variation, and the two drivers behind the eleven predictors

**Pooled within-range correlation matrix of the spring-mean anomalies** (standardised per range, 1583 range-years):

| | Td | SWE | snowfall | LH | SH | SWdown | LWdown | T |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Td | 1 | −.10 | −.04 | −.40 | .27 | −.44 | **.87** | **.65** |
| SWE | | 1 | **.53** | .44 | **.50** | −.03 | −.09 | −.40 |
| snowfall | | | 1 | .40 | .45 | −.33 | .10 | −.38 |
| LH | | | | 1 | .24 | −.12 | −.20 | −.50 |
| SH | | | | | 1 | −.38 | .30 | −.18 |
| SWdown | | | | | | 1 | **−.66** | .02 |
| LWdown | | | | | | | 1 | **.57** |

PCA of the eight: PC1 37 % (T +.49, Td +.49, LWdown +.46, LH −.37, SWE −.30: *warm, humid, cloudy-sky spring with the
ground melting out*), PC2 33 % (SH +.49, snowfall +.44, SWE +.36, LWdown +.34, SWdown −.44: *snowy, overcast spring with
the ground still snow-covered*), PC3 10 %. Two axes carry 70 % of the predictor variance.

**Per range, the correlation of each spring-mean anomaly with the spring-mean T anomaly** (median and interquartile range
over 160 ranges): Td 0.82 (0.48 to 0.95), LWdown 0.66 (0.40 to 0.82), LH −0.60 (−0.82 to −0.30), SWE −0.44 (−0.66 to
−0.21), snowfall −0.43 (−0.68 to −0.13), SH −0.17 (−0.59 to 0.17), SWdown 0.01 (−0.26 to 0.38). Within a typical range a
warm spring is a humid, longwave-rich, low-snow spring; it is *not* a sunnier spring.

**The antecedent snowpack is set by winter snowfall, not winter temperature** (pooled, standardised): March SWE vs
winter snowfall sum 0.61; vs winter T −0.15. So "March SWE" is a precipitation-side driver, largely independent of the
temperature-side driver in cold ranges, and increasingly confounded with it toward the subtropics (per-range r(spring T,
March SWE) median −0.06 at 60–90°, −0.36 at 50–60°, −0.44 at 35–50°, −0.51 at 0–35°).

---

## 6. How much variance, and by what: pooled models, commonality, and the honest per-range number

**How.** Predictors standardised within range (so a coefficient is "standard deviations of onset per standard deviation
of the predictor" and the pooled fit is the within-range, fixed-effects relationship); OLS through all range-years;
95 % CI from a bootstrap that resamples *water years* (all ranges' values for a drawn year move together, 300 draws).
The **driver set** is fixed a priori: spring-mean T, March SWE (`spring_month_1`, the antecedent pack), spring snowfall
sum, spring SWdown and LWdown, plus dewpoint, winter snowfall and winter T in the full model. The state variables
(spring-mean SWE, LH, SH) are fitted separately to show what they do.

| pooled model (1583 range-years, 160 ranges) | R² | standardised coefficients [year-block bootstrap 95 % CI] |
| --- | --- | --- |
| spring T | 0.22 | −0.47 [−0.54, −0.41] |
| March SWE | 0.20 | +0.45 [+0.39, +0.50] |
| spring snowfall | 0.18 | +0.42 [+0.37, +0.48] |
| spring SWdown | 0.00 | +0.03 [−0.03, +0.09] |
| T + March SWE | 0.32 | T −0.37 [−0.45, −0.30]; SWE +0.33 [+0.26, +0.38] |
| T + March SWE + spring snowfall | 0.35 | T −0.30; SWE +0.28; snowfall +0.20 (all CIs exclude 0) |
| + SWdown + LWdown | 0.36 | T −0.35; SWE +0.25; snowfall +0.24; SWdown +0.17 [+0.08, +0.26]; LWdown +0.09 [−0.04, +0.21] |
| all eight drivers | 0.38 | T −0.31; Td −0.17; SWE +0.17; snowfall +0.25; winter snowfall +0.14; SWdown +0.20; LWdown +0.24; winter T −0.03 [−0.12, +0.07] |
| *state:* spring-mean SWE | 0.27 | +0.52 |
| *state:* spring LH | 0.17 | +0.41 |

**Commonality analysis** of T + March SWE + spring snowfall (R² = 0.35): unique to T 0.07, unique to SWE 0.06, unique to
snowfall 0.03, **shared 0.18**. Half the explained variance cannot be assigned to one driver: warm springs are low-snow
springs.

Two things to note and not over-read. The all-driver model reaches 0.38, so ERA5-Land monthly anomalies explain about
40 % of the within-range interannual onset variance globally, and the first two drivers get you 0.32. And SWdown turns
weakly *positive* once T and snowfall are held fixed (+0.17): at the same temperature, a sunnier spring is a drier,
clearer-sky spring with lower longwave and less condensation (a suppression effect through the cloud axis), not a
melt-accelerating one. It is a curiosity at this stage, not a result.

**By latitude band**: the table in the summary (item 4). Poleward of 60°, T + March SWE gives R² = 0.41 with SWE at
+0.58 and T at −0.26; between 35° and 60°, T leads (−0.43 to −0.49) and SWE adds (+0.30); in the tropics and subtropics
snowfall leads (+0.32) and nothing reaches 0.25.

**The honest per-range number: leave-one-year-out predictive R²** (133 ranges with n ≥ 8; negative means the model
predicts a withheld year worse than the range's mean does):

| model | 25th pct | median | 75th pct | share of ranges improved over T alone (by > 0.1) |
| --- | --- | --- | --- | --- |
| spring T | −0.37 | **0.02** | 0.31 | |
| T + March SWE | −0.46 | 0.16 | 0.43 | 0.47 (0.35) |
| T + spring snowfall | −0.48 | 0.00 | 0.32 | 0.32 (0.23) |
| T + SWdown | −0.62 | −0.18 | 0.26 | 0.25 (0.14) |
| T + March SWE + snowfall | −0.74 | 0.08 | 0.40 | |

For the median range, the spring-temperature line does not predict a withheld year better than the range's mean; in
the 47 ranges with r² ≥ 0.5 the median predictive R² is 0.38. Adding a noise predictor (SWdown) makes things worse in
three ranges out of four, which is the yardstick: SWE beats that yardstick clearly, snowfall barely. Three-predictor
per-range fits blow up (min −95) on high-leverage years. **This is the statistical reason the framework in section 8
builds pooled and hierarchical models instead of per-range multivariate ones.**

```python
# scratch/climate_drivers/d_multivariate.py, the two pieces that matter
Dz[cols] = D.groupby(level=0)[cols].transform(lambda c: (c - c.mean()) / c.std(ddof=0))   # standardise within range
def year_block_bootstrap(df, xs, n=300):                                                  # resample WATER YEARS
    return np.percentile([ols(pd.concat([df[df.index.get_level_values(1) == y] for y in
                          rng.choice(years, len(years), replace=True)]), xs)[0] for _ in range(n)], [2.5, 97.5], axis=0)
def loo_r2(g, xs):                                                                        # leave-one-year-out
    ...; press += (y[i] - A[i] @ lstsq(A[~i], y[~i]))**2; return 1 - press / ((y - y.mean())**2).sum()
```

---

## 7. The weakly temperature-sensitive ranges, by regime

46 ranges with n ≥ 8 have spring-T r² < 0.2. Under the a priori driver set (fixed windows, no month search) they sort
into four groups (`scratch/climate_drivers/per_range_models.csv` has every range):

**(a) Polar and subpolar, snow-amount paced** (r² for March SWE, r positive): Kolyma 0.64, Anadyr 0.65, Kamchatka 0.66,
Arctic Ocean 0.57, East Siberian Plateaus 0.49, Koryak 0.48, Yukaghir 0.47, Chersky 0.41, Alaska Range 0.40, Northern
Scandes 0.41 (its T r² is 0.30), Verkhoyansk 0.40. The T + SWE + snowfall leave-one-year-out R² is 0.6–0.8 for Kamchatka,
Koryak and the Alaska Range: these are genuinely predictable ranges whose driver is not temperature. Physically: near the
solstice the energy for melt is abundant and reliable every year, so the year-to-year timing is set by how much snow must
be ripened (cold content and liquid-water-holding capacity both scale with SWE) and by how late the last snowfalls came.
The Brooks Range is the exception in this group: nothing explains it (SWE 0.13, n = 9).

**(b) Dry subtropical ranges where the amount of snow decides whether there is a melt season at all**: Southwest Basins
and Ranges (SWE 0.69, spring snowfall 0.72, n = 9), High Atlas (SWE 0.40), Gobi-Altai (0.43), Mongolian Altai (0.38),
Kazakh Uplands (0.44), Tasmania (0.36). Plausible and consistent, but n is small and these ranges have few pixels.

**(c) Nothing in ERA5-Land explains them** (no a priori driver reaches r² 0.35, and no sign is consistent): the tropical Andes (all three cordilleras), Hawaii,
the interior Chinese and Tibetan ranges (Qin Ling, Min, Tanggula, Bayan Har, Altyn-Tagh, Alashan, Hengduan, Qionglai),
New Zealand (North Island, Kaikoura, South Island Hills), Hokkaido, the Adirondacks, Jura, Massif Central, Central
Uplands, Western Mediterranean Islands, Central Labrador, Dzhugdzhur. Three sub-stories: (i) low, marginal-snow ranges
where the onset anomaly is a small-number statistic of few pixels; (ii) the tropical Andes, where "onset" is a
late-season artefact by the pipeline's own rule and the seasonal cycle is weak; (iii) the dry Tibetan-margin ranges,
where sublimation, wind and a thin, cold, patchy pack make the S-1 wet-snow signal noisy and ERA5-Land's 9 km
temperature and modelled SWE unrepresentative. In none of these is there evidence for a *different ERA5 variable*.

**(d) Positive temperature sensitivity**: Verkhoyansk (T r² 0.26 with a positive slope; SWE 0.40; LWdown 0.38). The partial
correlation of T given SWE is still +0.51, so the positive sign is not a snow-amount artefact; a warm spring at 64°N in
the world's coldest inhabited region may mean a snowier spring (more precipitable water) or a later, wetter ripening.
One range, eleven years: a discussion point, not a finding.

**Exemplars where T is strong, for contrast**:

| range | n | r² T | r² March SWE | r² spring snowfall | r(T, SWE) | partial r T\|SWE | partial r SWE\|T | LOO R² T | LOO R² T+SWE |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sierra Nevada | 11 | 0.71 | 0.63 | 0.64 | **−0.77** | −0.60 | +0.43 | 0.51 | 0.41 |
| Colorado Plateau | 11 | 0.64 | 0.70 | 0.40 | **−0.87** | −0.25 | +0.49 | 0.49 | 0.53 |
| Cascade Range | 11 | 0.75 | 0.50 | 0.30 | −0.34 | −0.94 | **+0.87** | 0.45 | **0.83** |
| Tian Shan | 11 | 0.41 | **0.84** | 0.04 | −0.71 | **+0.04** | +0.85 | 0.00 | **0.70** |
| Himalaya | 11 | 0.89 | 0.06 | 0.60 | −0.35 | −0.94 | −0.28 | 0.75 | 0.73 |
| European Alps | 11 | 0.57 | 0.05 | 0.35 | −0.04 | −0.77 | +0.30 | 0.35 | 0.19 |
| Qilian Mountains | 11 | 0.56 | 0.08 | 0.15 | −0.52 | −0.73 | −0.19 | 0.24 | 0.11 |

The Sierra Nevada and the Colorado Plateau are *both* temperature- and snow-amount-sensitive and the eleven years
cannot tell the two apart (partial correlations both moderate, r(T, SWE) below −0.75). The Cascades and the Tian Shan
can: in the Cascades both drivers act and the two-predictor model nearly doubles the predictive R²; in the Tian Shan the
whole temperature signal is explained by the snow amount (partial r of T given SWE is 0.04). The Himalaya and the Alps
are temperature-paced with little snow-amount effect. So "strong T sensitivity" in the map is a mixture of two
mechanisms, and the worked example should say so.

---

## 8. What is left out

**Processes the eight monthly anomalies cannot carry**

- *Threshold crossing, not monthly mean.* Ripening starts when the surface first stays near 0 °C for days; a monthly mean
  hides when in the month that happened. Evan & Eisenman (2021) make exactly this point: the sensitivity depends on where
  the annual temperature cycle sits relative to the melt threshold, which is why continental interiors (cycle crosses 0 °C
  steeply) are less sensitive than maritime ranges. Positive degree-days before onset, or the date of the first n-day
  run above 0 °C, is the natural predictor and needs daily data.
- *Absorbed, not downwelling, shortwave.* Snow reflects most of SWdown; the absorbed part is governed by albedo (grain
  size, age, dust and black carbon, canopy shading), none of which is in the store. ERA5-Land has `snow_albedo`,
  `forecast_albedo` and `surface_net_solar_radiation_sum` (modelled, but at least the right quantity).
- *Rain-on-snow and liquid precipitation.* `total_precipitation_sum − snowfall_sum` is rainfall; its timing (daily) is
  what matters.
- *Cold content and pack structure.* Snow temperature (`temperature_of_snow_layer`), density, depth; the energy needed
  to bring the pack to 0 °C scales with SWE × temperature deficit. This is the mechanism behind the polar SWE result and
  it is testable with the model's own snow temperature.
- *Wind* (turbulent exchange, sublimation, redistribution): `u/v_component_of_wind_10m`.
- *Cloud cover and its diurnal timing* (only through SW/LW proxies).
- *The sensing floor.* Onset is the Sentinel-1 backscatter minimum with a 6–12 day revisit; the per-pixel `temporal_resolution`
  in the pyramid is a hard floor on the detectable anomaly (several days) and should bound the explainable variance
  per range before any "unexplained" is attributed to physics.
- *The response itself.* Onset ≈ the start of wet-snow/runoff, not melt-out (Evan & Eisenman's ζ) nor peak SWE; the
  drivers of the three differ, and ERA5-Land's `snowmelt_sum` and `runoff_sum` are the model's own analogue of ours.
- *Upstream drivers* (ENSO, PDO, PNA, atmospheric rivers) act through T and snowfall; `teleconnections.ipynb` has the
  indices but they have not been regressed against the range-years.
- *Land cover and canopy*: FCF ≤ 50 % masks the densest forest; within the kept pixels canopy still shades and intercepts.

**ERA5-Land monthly bands that exist and are not acquired** (checked live today against `ECMWF/ERA5_LAND/MONTHLY_AGGR`,
150 bands): `snowmelt_sum`, `snow_cover`, `snow_albedo`, `snow_density`, `snow_depth`, `snow_evaporation_sum`,
`temperature_of_snow_layer`, `skin_temperature`, `temperature_2m_max` and `_min` (monthly extremes of the daily values,
a cheap proxy for threshold crossing), `total_precipitation_sum`, `surface_net_solar_radiation_sum`,
`surface_net_thermal_radiation_sum`, `forecast_albedo`, `u_/v_component_of_wind_10m`, `surface_pressure`,
`runoff_sum` / `surface_runoff_sum` / `sub_surface_runoff_sum`, soil temperature and moisture at four levels, and the
`_max/_min` of every accumulated variable. Adding bands is cheap in the current design (one variable × year × month
chunk each; the ledger handles refills); daily-derived indices need the daily collection (`ECMWF/ERA5_LAND/DAILY_AGGR`)
and a new per-water-year job that reduces days to a handful of indices per cell.

---

## 9. A framework for the methodical investigation

The principle: **hypotheses first, windows and variables fixed before looking, inference pooled across ranges under a
dependence-preserving null, per-range results descriptive.** Six phases; each ends in a results table with provenance
and, where noted, one figure. Phases 0–2 need nothing new from Azure or Earth Engine; 3 and 5 need a modest data
extension.

**Phase 0 — Foundations in the existing notebook (a day).**
- Declare each variable's role (driver / antecedent state / consequence) in a table in `temperature_sensitivity.ipynb`,
  and remove `total_surface_energy_sum` and `total_radiation_sum` from `HEATMAP_VARIABLES` (mis-specified, section 2).
- Add the *selection null* next to the heatmaps: the best-of-99 distribution under year-shuffling (section 3) as a one-line
  caption number, and a cross-range summary panel (median r² by variable × month, with the 1/(n−1) line) in front of the
  eleven per-range heatmaps. The heatmaps stay as the supplementary sweep.
- Compute the sensing floor per range from the pyramid's `temporal_resolution` (median revisit over the range's pixels)
  and store it in the metrics table; report "explainable variance" as 1 − floor²/var(onset anomaly).
- Write `results/<version>/climate_driver_screen.csv`: (range, variable, window, n, slope, r, r², p), and
  `climate_driver_pooled.csv`: the section-6 pooled models with year-block bootstrap CIs, by band. Port the code from
  `scratch/climate_drivers/`.

**Phase 1 — Pre-registered driver set and windows.** Drivers: spring-mean T; antecedent SWE (the month before the range's
median onset month, falling back to `spring_month_1`); spring snowfall; winter snowfall; winter T; LWdown; SWdown (kept
as the negative control). For the 30 early-melting ranges, the window is onset-centred (section 4). State variables
(spring SWE, SH, LH) are reported as diagnostics only. Every hypothesis below names its predictor, window, expected
sign, and test *before* the fit:

| hypothesis | predictor(s) | expected sign | test |
| --- | --- | --- | --- |
| H1 thermal pacing weakens poleward | spring T | − ; \|coef\| decreasing with \|lat\| | pooled coefficient by band, year-block bootstrap; Spearman of per-range slope vs \|lat\| |
| H2 snow-amount pacing strengthens poleward | antecedent SWE | + ; increasing with \|lat\| | as H1 |
| H3 late-season snowfall delays onset independently of T | spring snowfall given T | + | pooled partial coefficient; commonality |
| H4 winter warmth advances boreal onset (mid-winter melt, rain-on-snow) | winter T, winter Td, rainfall (Phase 3) | − in 50–60° only | band-specific pooled model |
| H5 downwelling shortwave anomalies do not pace onset | spring SWdown | 0 | negative control: CI includes 0 in every band |
| H6 T and SWE effects are separable where they are uncorrelated | r(T, SWE) per range | — | compare partial correlations in low- vs high-confounding ranges |
| H7 the sensitivity depends on the annual cycle's position relative to 0 °C (Evan & Eisenman) | T₀ − Tₘ, T₁ per range (the E&E comparison notebook already fits T₀/T₁ globally) | per-range T slope more negative where \|T₀ − Tₘ\| is small and T₁ small | regression of per-range slope on T₀, T₁ across ranges (160 points, no small-n problem) |

**Phase 2 — One hierarchical model instead of 160 regressions.** A mixed-effects model on the 1580 range-years:
`onset_anomaly ~ T + SWE + snowfall + (T + SWE + snowfall | range)`, predictors standardised within range, with
*range-level* covariates explaining the slopes (|latitude|, T₀ − Tₘ, T₁, median onset, mean elevation, mean SWE). This
gives shrunken per-range slopes with honest intervals (partial pooling is exactly what n = 11 needs), a test of the
regime hypotheses as slope-by-covariate interactions, and one cross-validation scheme: leave-one-*year*-out across the
whole model (never leave-one-range-out, the ranges are the population). `statsmodels.MixedLM` suffices for a first pass;
`bambi`/PyMC if posterior intervals are wanted. Output: `climate_driver_mixed_model.csv` and a two-panel figure (slope
of T and of SWE against |latitude|, shrunken estimates with intervals, coloured by regime).

**Phase 3 — Extend ERA5-Land by the bands that carry the missing processes.** Monthly first (cheap, same ledger):
`surface_net_solar_radiation_sum`, `surface_net_thermal_radiation_sum`, `total_precipitation_sum` (→ rainfall),
`temperature_2m_max`/`_min`, `snow_albedo`, `snow_cover`, `temperature_of_snow_layer`, `snow_density`,
`snowmelt_sum`, `u/v_component_of_wind_10m`. Then daily-derived indices per water year from the daily collection,
reduced in the acquisition job to one value per cell: positive degree-days in the 60 days before the cell's median
onset, date of the first 5-day run of T_max > 0 °C, rain-on-snow days (rain > 5 mm with SWE > 0), days since last
snowfall > 5 mm at the median onset, and the model's own melt-onset date (first day `snowmelt` > threshold). The
zonal-mean machinery already handles any variable on (water_year, month or index, lat, lon).

**Phase 4 — Within-range structure, the cleanest design available.** The cube already has the onset anomaly by
elevation × aspect × insolation class. Regress each *bin's* anomaly series on the *range-level* predictors: the
predictor is shared, the response varies, so differences in slope across bins are not confounded by predictor noise.
Expected: low-elevation and sunny bins T-paced, high-elevation and shaded bins SWE-paced within the same range; the
slope-vs-elevation profile per range ties to the lapse-rate analysis and to E&E21's threshold-crossing argument (the
bins whose spring T climatology straddles 0 °C should be the most sensitive). Needs the absolute temperature
climatology per bin (ERA5-Land acquisition store + bin elevations), not just anomalies.

**Phase 5 — Ground truth at stations.** The SNOTEL/CCSS daily archive is cached and credential-free. Per station and
year: peak SWE, the melt-onset date from SWE (first sustained decline after peak) and, where available, from soil
temperature/moisture; positive degree-days; rainfall. Repeat the H1–H4 regressions with *measured* drivers at the western
SNOTEL/CCSS stations, and compare with the ERA5-Land-driven results for the ranges that contain them. This calibrates the
one modelled driver (ERA5-Land SWE against measured peak SWE anomalies) and checks whether the S-1 onset behaves like the
station melt onset with respect to the drivers. Reuse `evan_eisenman_2021_comparison.ipynb`'s station loader.

**Phase 6 — Reporting.** One notebook, `climate_drivers.ipynb`, that reads the cube and the extended zonal means, writes
the three results tables above, and produces (i) the variable × month cross-range summary with its null line, (ii) the
band × driver pooled table, (iii) the hierarchical-model slope panels, (iv) a driver-regime world map (each range
coloured by its dominant shrunken driver: temperature / snow amount / snowfall / unresolved, with the sensing floor
as transparency). The current eleven heatmaps move to a supplementary sweep folder like the scatterplot sweeps.
Statistical rules carried over from the stats note: year-block bootstrap for every pooled interval, the shared-year
permutation for every "how many ranges" statement, Theil–Sen next to OLS, no per-range R² in the main text.

**Decisions to make first (yours).**
1. Whether the paper's claim stays "spring temperature sensitivity" with a two-mechanism caveat (T *and* snow amount,
   inseparable in the Sierra), or becomes "two regimes: thermally paced mid-latitudes, snow-amount paced high
   latitudes". The pooled numbers support the second; it changes the framing of the world map.
2. Whether to spend the Earth Engine budget on the monthly bands (a day, no new machinery) or go straight to daily
   indices (a new job type, the bigger scientific gain).
3. Whether Phase 4 (bin-level) or Phase 5 (stations) is the validation the reviewers will want; both are cheap, the
   station one is the more conventional.

---

## 10. Caveats specific to this note

- All ERA5-Land snow quantities are model output without land data assimilation; the SWE result is a statement about
  the *modelled* accumulation-season snowpack, which tracks winter snowfall (0.61) and is known to be biased in complex
  terrain. Phase 5 is the check.
- The n ≥ 6 rule for fits and the n ≥ 8 rule for the per-range model table drop the least-sampled ranges; the pooled
  models use all 1583 range-years with complete predictors.
- The onset-relative windows (section 4) use the range's median onset month; for ranges whose onset spans two months
  across elevation the single month is a compromise.
- The selection null shuffles onset years within each range independently, which breaks the shared-year structure;
  it is the right null for "is the best cell noise", not for field significance (that is the stats note's section 4).
- The band tables split at 35°, 50° and 60° of centroid latitude, as in the stats note; the 0–35° band mixes the
  tropical Andes with the Himalaya and Tibet.

---

## 11. Addendum (2026-09-08, later): the high-latitude mechanism, more variables, and the aspect / clear-sky test

Added after the first read-through. Scripts: `scratch/climate_drivers/e_aspect_pilot.py` (aspect amplitude and phase per
range and per range-year; `aspect_contrast_climatology.csv`, `aspect_contrast_by_year.csv`,
`aspect_spread_definitions.csv`).

### 11.1 What we think is happening poleward of 60°

**The observations (v10, 28 ranges).** They melt in April (13) or May (14), under 16–24 h of daylight. Their interannual
variability is not smaller than the mid-latitudes' (median MAD 11.5 days in both bands), so the energy supply is not
"reliable enough to remove the variance"; the variance is paced by something else:

| median r² (r), >60° ranges | Dec | Jan | Feb | Mar | Apr | May | MAM mean | onset month |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SWE | .34 (+.58) | .33 | .42 | .47 (+.69) | .48 | **.59 (+.77)** | .54 | .58 |
| T2m | .12 | .07 | .04 | .05 (+.01) | .26 (−.42) | **.29 (−.53)** | .22 | **.37** |
| LWdown | .08 | .05 | .02 | .05 | .09 | **.24 (−.44)** | .21 | .12 |
| SWdown | .02 | .07 | .11 | .04 | .08 | .09 (**+.26**) | .09 | .07 |
| snowfall | .10 | .09 | .04 | .15 | .08 | .15 (+.39) | .18 | .10 |

Three things stand out. The SWE anomaly predicts the onset from *December* onward (r² 0.34 rising to 0.59 in May): the
accumulation trajectory sets the date months ahead. The temperature anomaly matters only in the onset month itself
(March 0.05, April 0.26, May 0.29; the onset-month r² is 0.37 versus 0.22 for the fixed MAM mean, which dilutes it with
an irrelevant March). And in May, more downwelling longwave goes with earlier onset (r = −0.44) while more downwelling
shortwave goes with *later* onset (r = +0.26): a clear May sky at 65°N is a cold, dry air mass; a cloudy one is warm
advection. T and SWE are uncorrelated within these ranges (median −0.06), so their effects add: pooled R² = 0.41 with
SWE at +0.58 and T at −0.26 (section 6).

**The mechanism we propose.** *Snow-amount-limited ripening with a synoptic trigger.* (i) Through March the air over
these ranges is far below 0 °C, so a temperature anomaly does not move the surface across the melt threshold; this is
the threshold argument of Evan & Eisenman (2021) in interannual form, and it is why the spring-mean sensitivity is weak
there rather than because temperature is unimportant. (ii) The energy needed to bring the pack to runoff is its cold
content plus the liquid-water deficit, both proportional to SWE; cold content is built by snowfall, not by the energy
balance (Jennings, Kittel & Molotch 2018, The Cryosphere 12, 1595), and shallower snow melts earlier (Musselman et al.
2017, Nature Climate Change 7, 214). With near-continuous daylight the May energy supply is large and repeatable, so
*how much snow must be ripened* sets the date more than the supply anomaly does. (iii) Within the onset month the
trigger is a warm, moist, cloudy air-mass intrusion that raises downwelling longwave and 2 m temperature together, as
Mortin et al. (2016, GRL 43, 6636) showed for Arctic sea-ice melt onset and Mioduszewski et al. (2014, JGR Atmospheres)
for land snow in northern Canada; our May LWdown and SWdown signs are that pattern. (iv) Caveats: ERA5-Land SWE is a
model product (no land data assimilation) and tundra snow is wind-redistributed and sublimates; under polar day there is
little nightly refreeze, so the Sentinel-1 backscatter minimum may behave differently than at mid-latitudes.

The three regimes, in one line each: **>60°** snow-amount-limited, synoptically triggered; **35–60°** temperature-paced
with the snow amount as modifier (and, at 50–60°, winter warmth and rain-on-snow); **<35°** snowfall-event-paced and
mostly below the noise floor.

### 11.2 Variables to bring in, in priority order

| priority | variable / index | why | source (checked in Earth Engine 2026-09-08) |
| --- | --- | --- | --- |
| 1 | positive degree-days in the 30 days before the range's median onset; date of the first 5-day run with T_max > 0 °C | the threshold-crossing predictor the monthly mean hides | `ECMWF/ERA5_LAND/DAILY_AGGR` (`temperature_2m_max/min`) |
| 2 | peak SWE (max Dec→onset), snow depth, density, `temperature_of_snow_layer`; cold-content proxy = SWE × (0 − T_snow) | the high-latitude driver, and the E&E-style mechanism test | ERA5-Land monthly and daily |
| 3 | onset-month LWdown and dewpoint; count of "intrusion days" (daily LWdown anomaly > +1 σ) | the synoptic trigger | ERA5-Land daily |
| 4 | rainfall = `total_precipitation_sum − snowfall_sum`; rain-on-snow days (rain > 5 mm with SWE > 0) | the 50–60° winter-warmth pathway | ERA5-Land daily |
| 5 | cloud fraction in the melt window | the clear-sky hypothesis (11.3) and the polar trigger | `NOAA/CDR/PATMOSX/V53` (daily, 0.1° — the ERA5-Land grid, 1979–); `MODIS/061/MOD10A1` class 250 (cloud-obscured days at 500 m on the snow pixels themselves; Terra 10:30 vs `MYD10A1` Aqua 13:30 gives the diurnal asymmetry); `MOD08_M3` `Cloud_Fraction_Mean_Mean` (1°, coarse) |
| 6 | albedo → absorbed shortwave = (1 − α)·SWdown | the term the current SWdown misses | `MODIS/061/MCD43A3` `Albedo_BSA_shortwave`, `MOD10A1` `Snow_Albedo_Daily_Tile`; ERA5-Land `surface_net_solar_radiation_sum` |
| 7 | 10 m wind (speed, prevailing direction in the melt window) | sublimation, redistribution, the lee-vs-sun phase test | ERA5-Land |
| 8 | Arctic Oscillation (spring) next to Niño 3.4 / PDO / PNA | the upstream driver of the polar intrusions | NOAA PSL (as `teleconnections.ipynb`) |
| 9 | sensing covariates: per-pixel `temporal_resolution`, ascending/descending share | the noise floor and the refreeze-timing bias | the pyramid |

### 11.3 The aspect contrast, Himalaya vs Hindu Kush, and how to test "clear skies"

**First, the size of the effect to be explained.** The draft says the western-HMA ranges show > 40 days of aspect spread
and the Himalaya ~10. Recomputed on v10 under the analyses' view (bins ≤ 100 pixels masked, CHILI folded), with five
definitions of "spread":

| definition | Hindu Kush | Himalaya | ratio | all-range median |
| --- | --- | --- | --- | --- |
| max − min of the elevation-relative onset over all cells | 27.4 | 20.1 | 1.4 | 22.1 |
| the same after the triplets' all-aspects-present mask | 27.2 | 20.1 | 1.4 | 21.4 |
| max over elevation bands of the per-band aspect spread | 24.7 | 19.2 | 1.3 | 19.8 |
| median over elevation bands of the per-band spread | 12.8 | 9.0 | 1.4 | 10.7 |
| **count-weighted aspect curve, max − min** | **11.0** | **7.5** | **1.5** | 8.8 |
| cosine amplitude of that curve [days] | 5.4 | 3.6 | 1.5 | 3.9 |

The contrast is real but 1.3–1.5×, not 4×. The extreme-cell definitions are driven by thin bins: the Karakoram's 70-day
max − min comes from two cells of 180 and 318 pixels; the Transhimalaya's 25-day max − min is noise (count-weighted
amplitude 0.4 days, cosine fit r² 0.2). The number quoted in the draft should be re-derived from v10 with the
count-weighted metric before the hypothesis is tested against it.

**Second, a sharper observable: the phase.** Fitting onset(aspect) = a + A·cos(aspect − φ) to the count-weighted curve
gives the azimuth of earliest onset. Across 155 ranges it clusters at the equator-facing direction (median offset +6°,
IQR −16° to +27°), and the departures are informative:

| group | ranges (azimuth of earliest onset) |
| --- | --- |
| western HMA, near due south | Hindu Kush 182°, Tian Shan 179°, Pamir 185°, Karakoram 191° |
| monsoon-influenced, rotated toward east | Himalaya 163°, Tanggula 167°, Kunlun 159°, Hengduan Shan 138°, **Bayan Har 79°** |
| no aspect signal at all | Transhimalaya, Qiangtang (amplitude 0.3–0.4 days) |
| rotated toward west | Great Basin 231°, Northern Scandes 232°, Southern Scandes 219°, Iceland 209°, Meseta Patagónica 275° (SH; sun at 0°) |
| large-amplitude polar | Brooks Range 9.3 days at 186°, Mackenzie 8.8 at 191° (versus Northern Scandes 2.5, Iceland 2.6) |

An eastward rotation is what afternoon convective cloud predicts: if the sky is clear in the morning and overcast after
noon, the "sunny" aspect is east-southeast. A westward rotation is what afternoon warmth (or west-wind lee loading)
predicts. The Brooks–Mackenzie vs Scandes–Iceland pair is the polar version of the Hindu Kush–Himalaya pair: continental
clear Arctic spring vs maritime cloud, at a sun angle where the direct-beam aspect contrast is largest.

**Third, the tests, from cheapest to most decisive.**

1. *Cross-range regression (no small-n problem).* Response: the count-weighted amplitude per range. Predictors: the
   melt-window cloud fraction over the range's snow pixels (PATMOS-x daily 0.1° or MOD10A1 cloud-obscured fraction, in
   the range's onset month ± 1) and the *clear-sky direct-beam contrast* between the equator- and pole-facing 30° slopes
   at the range's latitude and onset date (a solar-geometry integral; it removes the latitude/sun-angle confound that
   makes the Alps' 9-day amplitude unremarkable). Prediction: amplitude ≈ k · (1 − cloud fraction) · clear-sky contrast.
   The Hindu Kush–Himalaya pair, and the Brooks–Scandes pair, should sit on one global line.
2. *Direct-beam fraction instead of cloud fraction.* ERA5 single-level `fdir`/`ssrd` (CDS, not in Earth Engine) gives the
   direct fraction per cell and hour; the aspect amplitude should scale with the melt-window direct fraction. This is the
   physically exact predictor of the hypothesis and worth one CDS download for the ~160 range polygons.
3. *The phase test (sun vs wind).* Compare φ per range with the solar-noon azimuth shifted by the afternoon-heating lag
   and with the melt-season prevailing wind direction + 180° (lee). Radiative control keeps φ near south whatever the
   wind; depositional control tracks the lee. For the eastward-rotated monsoon ranges, the Terra (10:30) vs Aqua (13:30)
   cloud-obscured fraction difference from MOD10A1/MYD10A1 should be largest where φ rotates furthest east.
4. *The within-range interannual test — piloted here, null so far.* The year's aspect amplitude (from the cube's per-year
   aspect axis) against the year's spring SWdown anomaly, pooled over 1396 range-years standardised within range:
   coefficient +0.07 [−0.01, +0.13], R² 0.005; by band 0.00–0.01; per range in HMA, r(amplitude, SWdown) between −0.21
   (Hindu Kush) and +0.72 (Pamir). With the monthly SWdown anomaly as the cloud proxy and noisy per-year amplitudes, the
   design has no power; redo it with the onset-month cloud fraction before drawing a conclusion. If it stays null, the
   cross-range test carries the hypothesis.
5. *The along-arc gradient.* Slice the Himalaya into longitude segments (the case-study path that reads the store, or the
   HydroBASINS level-6 units that fall inside the polygon) and plot amplitude and phase against the melt-window cloud
   fraction from west to east; two polygons become a monotonic gradient or do not.
6. *Confounders to rule out.* Solar geometry (35° vs 30° N; onset date), the elevation distribution of the snow pixels
   (Himalayan snow sits higher), glacier and debris-covered pixels, snow-depth differences (deep high-elevation snow
   melting inside the monsoon), and the elevation-dependent rotation the draft already notes (fit φ per elevation band).

---

## 12. Addendum: why the triplet panels look more different than the metric says, and how much the topographic controls move between years

Scripts: `scratch/climate_drivers/f_interannual_topo.py` → `interannual_topography.csv` (per range and water year: lapse
rate, its SE, aspect amplitude, phase), `interannual_topography_summary.csv`.

### 12.1 Hindu Kush vs Himalaya: what the panel shows and what the number measures

The v10 triplets (`figures/v10/triplets/`) draw the elevation-relative median onset on a ±15-day scale; the archived
v9 notebook drew single-range panels on ±40 and the composite grid on ±10. The per-elevation cosine fits explain the
visual contrast:

| | Hindu Kush | Himalaya |
| --- | --- | --- |
| band of strongest aspect control | 3900–4200 m, amplitude **12.3 d** (spread 24.7 d), phase 179° | 4100–4400 m, amplitude **5.5 d** (spread 11.5 d), phase 165° |
| share of the range's pixels in that band | ~2 % per 100 m bin | ~4 % per 100 m bin |
| pixel-weighted mean amplitude over all bands | 6.0 d | 3.8 d |
| amplitude of the pixel-weighted curve (phase-coherent) | 5.4 d | 3.6 d |
| phase coherence across elevation bands | 0.66 | 0.91 |
| sign reversal (north-facing earliest) | below ~2500 m and above ~5900 m | below ~2800 m and above ~6200 m |

So three things are going on. (i) **The eye reads the peak band and its colour saturation**: ±12 days fills 80 % of a
±15 scale, ±5.5 days a third of it, and the peak-band ratio is 2.2×, which is what the panels show. (ii) **The
whole-range number is pixel-weighted**, and the Hindu Kush's strongest band is pixel-poor (its pixel mass sits at
2500–3500 m where the amplitude is only 0.4–8 d and still rising), while the Himalaya's moderate band is pixel-rich; that
narrows the ratio to 1.5–1.6×. (iii) **Both ranges reverse sign at their low and high ends**; the Hindu Kush's reversals
cost more coherence (0.66 vs 0.91). The polar projection adds a fourth effect: area grows with radius, so the outer
(low-elevation) rings, where both patterns are weak or reversed, get more visual area than their pixels justify.

The draft's "> 40 days vs ~10 days" is reproducible only without a pixel threshold. On v10, max − min of the
elevation-relative onset over all cells: no threshold 45 d (Hindu Kush) vs **58 d (Himalaya)**, from 1–5-pixel cells;
CHILI classes kept and bins ≤ 100 px masked 38 vs 26; the analyses' view 27 vs 20; bins ≤ 1000 px masked 27 vs 17.5.
The Karakoram goes 207 → 70 → 20 d as the threshold rises: its "> 30 days" is thin-bin noise. Every thresholded
definition gives 1.3–1.6×; the honest statement is "about twice as strong in the band where each is strongest, about
1.5× over the whole range". The Himalaya's "rotation with elevation" is a drift from 165° (3000–4500 m) to 145–150°
(4700–5100 m, i.e. south-east) and back to 180° (5300–5900 m), plus the end reversals both ranges share; the Hindu Kush
holds 177–183° from 2650 to 5700 m. The south-east excursion is the part that the afternoon-cloud hypothesis speaks to.
The low-elevation reversal (north-facing earlier than south-facing below ~2500–2800 m, few pixels) is unexplained and
worth a look on its own before it is called "rotation".

### 12.2 Interannual lapse rate, aspect amplitude and phase (103 ranges with ≥ 9 years)

Figures: `scratch/climate_drivers/figures/interannual_topographic_controls.png` (per-range strips: lapse rate, aspect amplitude, aspect
phase, one dot per water year in the repo's water-year colours, black diamond = the eleven-year median fields, ranges grouped
by continent north to south) and `interannual_aspect_roses.png` (one compass rose per range: arrow direction = azimuth of
earliest onset, length = aspect amplitude, one arrow per water year, black = median fields). Produced by
`g_plot_interannual_topography.py`; to be promoted into a notebook if kept.

Per range and water year: count-weighted OLS of the aspect-mean onset on elevation (days per 100 m, with its SE), and
the cosine fit of the count-weighted aspect curve (amplitude, azimuth of earliest onset). The yearly values use only
that year's valid bins (the 30 % rule), so they are noisier than the median-based ones; the fit SE is the yardstick.

| | median across ranges | interannual SD (median) | max − min across years (median) | single-year fit SE (median) | SD / SE |
| --- | --- | --- | --- | --- | --- |
| lapse rate [d / 100 m] | 4.0 | **1.0** (25 %) | 3.4 | 0.28 | **3.8** |
| aspect amplitude [d] | 4.3 | **1.6** (37 %) | 4.9 | 0.25 | **6.4** |
| aspect phase [°] | equator-facing | circular SD **19°** | | | 53 ranges < 20°, 13 ranges > 60° (no coherent signal) |

| range | lapse rate: median-year, min–max over years | amplitude: median-year, min–max | phase: median-year, circular SD |
| --- | --- | --- | --- |
| Sierra Nevada | 2.9; 1.9–5.0 | 5.7; 4.0–8.2 | 182°; 7° |
| Cascade Range | 4.5; 3.1–6.2 | 7.1; 1.9–10.0 | 186°; 20° |
| European Alps | 4.8; 3.8–6.0 | 9.1; 5.8–12.2 | 185°; 4° |
| Pyrenees | 5.6; 3.1–7.5 | 10.0; 6.1–14.1 | 189°; 6° |
| Hindu Kush | 4.2; 3.4–5.1 | 5.4; 4.5–7.6 | 182°; 3° |
| Himalaya | 3.0; 2.6–3.8 | 3.6; 1.8–5.9 | 163°; 9° |
| Tian Shan | 3.2; 2.7–3.6 | 4.7; 2.9–5.4 | 179°; 4° |
| Karakoram | 4.7; 4.1–5.2 | 4.4; 1.7–7.0 | 191°; 16° |
| Brooks Range | 3.0; 2.3–3.8 | 9.3; 6.3–10.5 | 186°; 7° |
| Kolyma Mountains | 1.8; 0.6–3.0 | 6.0; 2.2–6.8 | 177°; 6° |

**What it says.** The topographic controls are not fixed properties of a range. The lapse rate moves by a quarter and
the aspect amplitude by more than a third from year to year, four to six times the fit uncertainty, so the change is
real; the *direction* of the aspect control is the stable part (within ±20° in half the ranges). The 11-year median
panel is the centre of a distribution, and the per-range "lapse rate" in the metrics table should carry its
interannual range.

**What drives the year-to-year change.** Not the year's mean anomaly, to first order. Pooled within range
(standardised): r(lapse rate, onset anomaly) = −0.15 (earlier years slightly steeper; −0.23 poleward of 60°),
r(lapse rate, spring T) = 0.02, r(amplitude, onset anomaly) = −0.04, r(amplitude, spring T) = 0.01, r(phase deviation,
onset anomaly) = −0.12. The Sierra Nevada illustrates it: WY2015 (−36 d) has a lapse rate of 3.8, WY2016 (−10 d) 5.0,
WY2021 and WY2022 (−13, −10 d) 1.9. The one structural signal is r(lapse rate, amplitude) = **+0.24**: years that
differentiate the onset by elevation also differentiate it by aspect. A common "differentiation" factor is the natural
reading, and its two candidates are the direct-beam fraction (clear skies sharpen both the elevation and the aspect
contrast) and a thin snowpack (energy-limited melt-out exposes topography more than a deep pack does). That is a third
observable for the clear-sky test in 11.3: in the cross-range regression, use the *year* × range panel and ask whether the
cloud fraction explains the joint lapse–amplitude anomaly. It also means the 11 years per range are not wasted on the
aspect question even though the monthly SWdown proxy gave a null (11.3, test 4).

---

## 13. Addendum (2026-09-09): within a range, does a year's lapse rate go with its aspect control?

Scripts: `scratch/climate_drivers/h_lapse_vs_aspect.py` (→ `lapse_vs_aspect_per_range.csv`, `lapse_vs_aspect_range_years.csv`),
`i_plot_lapse_vs_aspect.py` (→ `figures/lapse_rate_vs_aspect_amplitude_by_year.png`). The scratch parquet inputs were
rebuilt into `scratch/climate_drivers/` with `build_pairs.py` (the session scratchpad had been cleared).

**How.** For each of the 103 ranges with ≥ 9 water years, the Pearson (and Spearman) correlation across its years
between the year's lapse rate (d / 100 m) and the year's aspect amplitude (d), both from section 12. Null: the
amplitude years shuffled within each range, 2000 times. Artefact check: the year's coverage (share of the range's
median pixels with a valid onset that year) against both. Pooled within-range regressions standardised per range with
year-block bootstrap intervals.

**Yes: a modest, widespread, positive relationship.**

| statistic | observed | shuffle null, 95th percentile | permutation p |
| --- | --- | --- | --- |
| median per-range r | **+0.26** (Spearman +0.29) | +0.07 | < 0.0005 |
| share of ranges with r > 0 | **0.74** | 0.58 | < 0.0005 |
| ranges with p < 0.05 | **20** (18 positive, 2 negative) | 9 | < 0.0005 |
| pooled standardised slope, amplitude on lapse rate | **+0.23** [+0.13, +0.32] | | R² 0.05 |

It is not coverage: per-range r(lapse rate, coverage) has a median of −0.08 and r(amplitude, coverage) +0.07, and the
pooled slope is unchanged with coverage in the model (+0.24). It holds in every latitude band (median r 0.18–0.30) and
whether the range's aspect control is weak, moderate or strong (0.25, 0.24, 0.34). The phase joins in: years whose
azimuth of earliest onset wanders from the range's usual direction are years with a weaker lapse rate (median r with the
absolute phase deviation −0.20) and a weaker amplitude (−0.24). So there is a "topographically differentiated year"
versus a "flat year" within a range, and the three controls move together.

**Exemplars.** Hindu Kush r = +0.84, Baetic System +0.90, High Atlas +0.83, South European Highlands +0.79, Klamath
+0.72, Tian Shan +0.59, Sierra Nevada +0.57, Cascades +0.51, Pyrenees +0.49, Alps +0.31; Himalaya +0.14 and Karakoram
+0.23 (weak); Arakan −0.81 and Honshu −0.63 are the two significant negatives (monsoon-margin and maritime ranges with
small amplitudes).

**What the joint anomaly tracks: little, in the monthly ERA5-Land anomalies.** Pooled standardised coefficients
(year-block bootstrap): lapse rate on spring SWdown −0.22 [−0.29, −0.13], on March SWE −0.13, on LWdown +0.12, on the
range's onset anomaly −0.15; amplitude on SWdown +0.09, on spring snowfall −0.10, on LWdown −0.09; everything at R² ≤
0.05. Two readings. The sunnier-spring sign is *opposite* for the two controls (shallower lapse rate, slightly larger
amplitude), so "clear skies" as a single monthly-mean factor does not produce the co-variation; a sunnier spring plausibly
advances the insolation-limited high elevations more than the temperature-limited low ones (compressing the gradient)
while sharpening the aspect contrast. And deeper March SWE flattens the gradient (−0.13). The common factor behind the
+0.26 is therefore either sub-monthly (the weather during each range's own melt window, which the fixed spring mean
does not see) or in the snow distribution (a thin, energy-limited year exposes topography in both dimensions) rather than
in the seasonal-mean forcing. The onset-window daily indices and cloud fraction of section 11.2 are the data to test it
with; the joint lapse–amplitude anomaly per range-year is the response variable.

---

## 14. Synthesis (2026-09-09): synoptic vs gradual springs, the literature, Evan & Eisenman tested, and the aspect-by-elevation view

Scripts: `scratch/climate_drivers/j_evan_eisenman_test.py` (→ `evan_eisenman_cross_range.csv`),
`k_aspect_by_elevation_profile.py` (→ `figures/aspect_by_elevation_profiles.png`), `variance_vs_sensitivity.csv` (from an
inline check). Literature read from `../global_snowmelt_runoff_onset_analysis_paper/literature/` (abstracts extracted
with pypdf; the folder's `paper/notes_literature.md` already argues for the radiative-vs-advective framing).

### 14.1 The interpretation, and its precedent in Lundquist et al.

The section-13 result (a year's lapse rate and aspect amplitude co-vary within a range, r ≈ +0.26, and off-axis phase
years are flat years) reads as: **when a spring is set by a synoptic event, melt begins nearly everywhere at once and
neither elevation nor aspect gets to act; when the spring is quiet and melt progresses gradually, insolation paces it
and both the elevation gradient and the aspect contrast are at their largest.** That is Lundquist, Cayan & Dettinger
(2004, *J. Hydrometeorol.* 5, 327) in interannual, range-scale, global form: in their Sierra Nevada record, "synchronous
springs" (a late-winter trough followed by a strong ridge and a > 12 °C rise in five days; 8 of 87 years) start streamflow
at all elevations and on north- and south-facing cirques together, while in the most elevation-dependent years
"periods of little or no storm activity … precede the onset of spring snowmelt, allowing elevation and aspect to
influence snowmelt as spring arrives gradually". Our data extend this from one range and streamflow to 103 ranges
and 80 m onset, with the lapse rate and the aspect amplitude as the two measured axes of "elevation and aspect
influence".

Lundquist & Flint (2006, *J. Hydrometeorol.* 7, 1199) is a different mechanism and should be kept distinct: the March
2004 event melted before the equinox, at low sun, so *shading* mattered more than in a late-March event of equal warmth
(2002). That is a sun-angle-at-melt-time effect, the one our seasonal arc (the continents' sunny–shaded difference vs
median onset) documents climatologically. Interannually the two mechanisms pull in opposite directions in an early
year: lower sun (more aspect contrast) versus synoptic forcing (less). The data show which wins where: pooled over all
ranges the year's amplitude is unrelated to the year's onset anomaly (r = −0.04), but in the Sierra Nevada it is −0.63
(earlier years, larger aspect contrast: Lundquist & Flint in their own range) and in the Hindu Kush −0.44, while in the
Alps it is −0.02 and in the Pamir +0.28. Both citations belong in the paper, attached to different sentences.

Also relevant, from the folder: Rhoades et al. (2026, *Sci. Adv.*, "snow-eater heat waves") characterise the synoptic
trigger itself for the western US (3–5-day events, three to five per melt season, doubling melt rates, increasing in area
and frequency and shifting earlier since the 1850s) — the natural citation for what a "flat year" is made of, and for
the warming-implication sentence (more flat years, not only earlier ones). Pflug & Lundquist (2020) show that the
*spatial pattern* of snow is repeatable between years while its amplitude is not, which is the snow-side analogue of our
finding that the aspect *direction* is stable (circular SD 19°) while the aspect *amplitude* varies by a third.

### 14.2 What the draft says now, and what to change

The draft (`paper/drafts/draft_4_yet_another_version.md`) already has the idea qualitatively: "These are cases in which
large-scale atmospheric forcing overwhelms typical local controls and produces near-synchronous shifts across elevation
and aspect" (citing Lundquist 2004) and, in Results, "synoptic-scale forcing produces spatially coherent, regional
events that can overwhelm typical local controls". What it lacks is the quantification by water year, and it carries
several statements the v10 numbers do not support. In priority order:

1. **Add the interannual result as a finding, not an aside.** One paragraph in Results: lapse rate SD 1.0 d/100 m
   (25 %) and aspect amplitude SD 1.6 d (37 %) between years, four to six times the fit error; direction stable (19°);
   lapse rate and amplitude co-vary (median r +0.26, 74 % positive, 20 of 103 significant vs 9 by chance; shuffle
   p < 0.0005). Then the Lundquist 2004 sentence, now as a tested generalisation. Candidate figures: the roses grid, the
   strip panels, or the lapse-vs-amplitude small multiples (section 12/13 figures).
2. **Re-derive every aspect number with a count-weighted metric** (sections 11.3, 12.1). "> 40 days" (western HMA) and
   "~10 days" (Himalaya) become "peak-band amplitude 12 vs 5.5 days (spread 25 vs 12), whole-range 1.5×". Table 1's "Max
   aspect diff." column is the thin-bin maximum and must be recomputed (the Karakoram's "> 30" is two ~200-pixel cells);
   the Methods sentence "maximum spread … across 8 aspect bins within each elevation band" becomes the cosine amplitude
   on 24 bins with ≥ 16 aspects and ≥ 100 (better 1000) pixels per bin. The abstract sentence "40 days … less than 10
   days" goes with it.
3. **"Rotate with elevation" → "reverse at both ends, drift in the middle".** Both the Hindu Kush and the Himalaya put
   north-facing slopes first below ~2500–2800 m and above ~6000 m; the Himalaya additionally drifts south-east at
   4500–5200 m. Say that, and flag the low-elevation reversal as unexplained (few pixels; possibly ephemeral snow on the
   sunny side leaving only shaded persistent snow, or a late re-accumulation signal).
4. **Two regimes for the temperature sensitivity, with the confounding stated.** Mid-latitudes temperature-paced with
   the snow amount as co-driver (Sierra Nevada r(T, March SWE) = −0.77: the two are inseparable there); poleward of 60°
   snow-amount paced with a synoptic longwave trigger (section 11.1). Replace "72 % of ranges |r| > 0.4" with the
   permutation statement from the stats note; add the Evan & Eisenman cross-range test (14.4) as the mechanism sentence.
5. **High latitudes are not "low variability".** Median anomaly SD 7.8 d and MAD 12.0 d poleward of 60° vs 9.5 / 11.6 d
   at 35–50°. "Polar … mountain ranges, which melt under strong, consistent solar inputs, show MAD values under 15
   days" should become "melt under a reliable energy supply, so their year-to-year timing follows the snow amount rather
   than the spring temperature" (Mioduszewski et al. 2015: insolation and longwave drive melt where onset is in May–June,
   advection where it is in March–April).
6. **Radiation anomalies do not pace the interannual anomaly.** Downwelling shortwave explains 0.00 of the pooled
   within-range variance in every band; insolation sets the *climatological* pattern (aspect, latitude, the seasonal arc),
   the *anomalies* are paced by temperature/humidity/longwave and by snow amount. This sharpens the "radiative vs
   advective" language `notes_literature.md` already recommends over "temperature-dominated vs radiation-dominated".
7. **Tropical Andes and interior Tibet: say "unresolved", not "sublimation".** Nothing in ERA5-Land explains them, their
   anomaly SD is the largest anywhere (30–42 d), and Fyffe et al. (2025) describe the Peruvian pack below 5000 m as thin
   and ephemeral with rapid snowfall–melt cycles: the S-1 minimum there is not the same event as a mid-latitude ripening.
8. **Table 1**: add the interannual range (min–max over years) next to the lapse rate and the aspect amplitude, and the
   Theil–Sen slope next to the OLS sensitivity.
9. **Warming-implication paragraph**: the draft's "migrate along the continuum toward more spatially heterogeneous,
   temperature-sensitive melt" needs the counter-weight that earlier melt is often *synoptic* melt (Rhoades et al.
   2026; our flat years), which is *less* heterogeneous; and Zheng et al. (2022) found the melt-onset advance
   decelerating since the 1990s and the Southern Hemisphere trending the other way. Say "two opposing tendencies" rather
   than one direction.

### 14.3 Other literature to connect

| paper (in the folder unless marked) | connects to |
| --- | --- |
| Mioduszewski et al. 2015 | our latitude regimes: advection-paced where onset is Mar–Apr, insolation + longwave where May–Jun; their PCA is the model for a cross-range PCA of our five metrics (recommended in `notes_literature.md`) |
| Zheng et al. 2022 (radiometer melt onset, 40 yr) | AO/NAO control of pan-Arctic onset → add the spring AO to `teleconnections.ipynb`; the hiatus and SH divergence caveats for the warming paragraph; their Feb-2018 SSW vignette parallels our WY2015/2022/2023 events |
| Cortés, Vargas & McPhee 2011 (Andes streamflow timing) | central-Chile timing correlates with precipitation/ENSO (r > 0.8) more than temperature: matches our subtropical band, where spring snowfall and March SWE lead and T is weak (Dry Andes: T r² 0.48 but snowfall 0.53) |
| Fyffe et al. 2025 (Peruvian Andes) | why the tropical Andes are unresolved (thin, ephemeral snow below 5000 m) |
| Rhoades et al. 2026 (snow-eater heat waves) | the synoptic trigger; the flat-year mechanism; trend context |
| Pflug & Lundquist 2020 (snow-pattern repeatability) | stable aspect direction, variable amplitude |
| Bayle et al. 2025; the Alps–Pyrenees SMOD record | complementary melt-out-date products for validating onset-vs-melt-out differences and for the Alps' stability |
| *not in the folder:* Mortin et al. 2016 (GRL), Mioduszewski et al. 2014 (JGR) | the longwave / warm-moist-intrusion trigger at high latitude (section 11.1) |
| *not in the folder:* Jennings, Kittel & Molotch 2018 (TC); Musselman et al. 2017 (NCC) | cold content set by snowfall; shallower snow melts earlier — the SWE mechanism |
| *not in the folder:* Harpold & Brooks 2018 (PNAS, humidity determines ablation) | the dewpoint/longwave axis (PC1 of our predictors) as the physical reason T, Td and LWdown act together |

### 14.4 Evan & Eisenman, tested across our ranges

**Their argument, in our terms.** The date the temperature cycle crosses the melt threshold shifts by ΔT divided by the
seasonal warming rate at the crossing, ω·√(T₁² − (T₀ − Tₘ)²). Where the annual cycle *grazes* the threshold (maritime,
small T₁, or T₀ close to Tₘ) the rate is small and a degree moves the date by 10–20 days; where the cycle crosses
steeply (continental interiors, large T₁, T₀ far from Tₘ) a degree moves it by 2–4 days. Their Eq. 2 adds the factor
(1 + Rₐ/Rₘ) for the disappearance date. So limited spring-T sensitivity is *predicted* for continental Siberia, interior
Asia and the Arctic, not because temperature is unimportant but because the seasonal cycle there does the crossing
almost regardless of the anomaly, and (our section 11.1) because in March the surface is too far below 0 °C for an
anomaly to reach the threshold at all.

**The test.** T₀ and T₁ from the E&E comparison notebook's cached ERA5 climatology (0.7°, sinusoid fit, cell-mean
elevation), sampled at range centroids; Tₘ = 0.18 °C and Rₐ/Rₘ = 0.34 as in the notebook; observed = the metrics
table's OLS spring-T slope (n ≥ 8).

| | value |
| --- | --- |
| ranges with a defined crossing at the cell | 84 of 134 |
| Spearman ρ, E&E-predicted dζ/dT₀ vs observed slope | **+0.51** (p < 0.001); vs observed r +0.48 |
| median observed slope / median E&E prediction / median first-order crossing shift | −5.0 / −7.0 / **−5.2** d per °C |
| tercile of the warming rate at the crossing: grazing → middle → steep | observed median slope −6.7 → −5.8 → **−1.5** d/°C; median r −0.72 → −0.71 → **−0.27**; share r < −0.6: 0.79 → 0.68 → 0.18 |
| within 30–55° only (56 ranges) | ρ = +0.37 (p = 0.005): not just the polar contrast |

The first-order crossing shift reproduces the observed median almost exactly, and the steep-crossing tercile is where
the insensitive ranges are (Verkhoyansk T₁ = 27 °C, Kolyma 21, Brooks 18; Sierra Nevada and Cascades T₁ ≈ 9 °C, predicted
−13.6/−13.8, observed −12.6/−10.6). **Caveat that must be fixed before quoting:** 50 ranges are "undefined" because the
0.7° *cell-mean* climatology never reaches 0 °C — 49 of them on the warm side (Mediterranean, Iberian, Iranian and
Californian coast ranges, New Zealand, the tropical Andes), i.e. the cell is too low and warm to represent the range's
snow. Those ranges have a median observed slope of −5.9 d/°C and r −0.61, so they are not physically insensitive; the
test needs T₀ adjusted to the range's mean snow-pixel elevation (a lapse-rate correction, or ERA5-Land at 9 km).

**Beyond E&E: other controls in the literature.** E&E themselves regress ζ on both T₀ and P₀ (precipitation) and make
Rₐ/Rₘ (accumulation vs melt rate) part of the sensitivity. Cortés et al. (2011) find precipitation and ENSO, not
temperature, control Andean timing; Mioduszewski et al. (2015) find advected energy where onset is early and
insolation + longwave where it is late; Zheng et al. (2022) tie regional onset to AO/NAO and ENSO/SAM; Rhoades et al.
(2026) to multi-day heat waves; Mortin et al. (2016) to moisture intrusions; Jennings et al. (2018) to snowfall-built
cold content; Harpold & Brooks (2018) to humidity. All of these are the "other factors" of the original question and
all map onto the two axes we found (warm-humid-longwave; snow amount and snowfall).

**Is it just low variance?** No. Weak-T ranges (r² < 0.2) have a median anomaly SD of 10.3 days against 8.9 for
strong-T ranges (r² > 0.5) and the same MAD (12.2 vs 11.1); the polar band has the *smallest* SD (7.8 d) but still
comparable to the mid-latitudes (9.5). Only 6 of 133 ranges have an interannual SD below 5 days (Hokkaido, Caucasus,
Northern Altai, Changbai, Tian Shan, Yan), and two of those still show r² ≈ 0.56. The weak ranges with the most variance
(tropical Andes 30–42 d, Qin Ling 21, Hawaii 18, Min 17, Bayan Har 16) are unexplained by anything, not quiet. The
sentence for the paper: "insensitivity is a property of the forcing regime and, in a minority of ranges, of the signal,
not of a lack of variability". (`variance_vs_sensitivity.csv`.)

### 14.5 Aspect control as an elevation profile: the proposed per-range view

`figures/aspect_by_elevation_profiles.png` prototypes it for the Hindu Kush, Himalaya, Sierra Nevada and Alps. Three
panels per range, elevation on the shared y-axis: (left) the share of the range's snow pixels per 100 m band, so the
reader sees where the mass is; (middle) the cosine amplitude of the aspect curve per band, one dot per water year in the
water-year colours, the eleven-year median fields as a line; (right) the azimuth of the earliest-melting aspect per
band as deviation from equator-facing, marker size scaling with that band-year's amplitude, so weak bands fade.

What it shows that the polar triplet cannot: the Hindu Kush is one broad coherent band (3500–5200 m, 8–12 d, due south,
tight year to year) with reversals below 2500 m and above 5900 m; the Himalaya has two weaker maxima (4000–4300 m at
5.5 d; 6800 m at 8 d but almost no pixels), a south-east drift at 4500–5200 m, and the same end reversals; the Sierra
Nevada's amplitude *grows* with elevation (3 → 8 d) and its phase rotates clockwise toward south-west/west with
elevation (afternoon heating at high elevation); the Alps peak at 1500–2500 m (11 d) and halve above 3000 m, with a
reversal only at the glacierised top. In every range the year-to-year spread is smallest where the pixels are and the
direction holds from year to year in the pixel-rich bands. This is the honest replacement for a single "aspect
difference" number: report the peak-band amplitude with its elevation, the pixel-weighted amplitude, and the phase
coherence; show the profile for the exemplar ranges and sweep it for all ranges in the supplement (as the triplets are
now). A compact two-panel version (amplitude and phase profiles, pixel share as line thickness) would fit the composite
map's inset budget.

### 14.6 The (range, water year) view: the distribution of aspect timings at each elevation

`scratch/climate_drivers/l_aspect_timing_by_elevation_year.py` → `figures/aspect_timing_by_elevation_<range>_WY<year>.png`
(one range-year) and `figures/aspect_timing_by_elevation_<range>_all_years.png` (the median fields and the eleven years).
Rendered for the Sierra Nevada (WY2015, WY2023, all years) and the Hindu Kush (WY2020, all years).

**The form.** The polar triplet unrolled: rows are 100 m elevation bands, columns the 24 aspect bins (N → E → S → W → N),
colour the onset minus the band's median across aspects (RdBu ±15 d, red = earlier; the triplets' convention), with a
black tick per band at the cosine-fit azimuth of earliest onset. For one range-year the figure adds the eleven-year
median fields beside it, the difference (year minus usual pattern, PuOr ±10 d), the year's pixel share per band, and
a fifth panel with the *absolute* onset (DOWY) against aspect as one line per band coloured by elevation, so the lapse
rate (line spacing) and the aspect curve (line shape) are read together. Bands with fewer than 16 aspect bins of
≥ 100 pixels, or under 30 % of the band's usual pixels that year, are blank.

**What it shows.** Sierra Nevada WY2015 (36 days early): the usual south-first pattern holds below 3000 m but the
strongest contrast moved *up* to 3000–4000 m and the earliest azimuth swung to south-west/west there (+8 to +12 d on
north/north-east faces, −8 d on west faces relative to the usual year), while the lower bands lost aspect contrast: a
warm, dry, early year in which low elevations melted out together and the high country was where insolation still
differentiated. The Hindu Kush grid shows the same south-facing band 3000–5500 m in every year, with the year-to-year
change in its *intensity* (WY2020, WY2023 strong; WY2017, WY2024 pale), a thin band above 5800 m that flips between
years, and WY2018 with the low elevations missing (the 30 % coverage rule). The absolute-onset lines make the lapse
rate and the aspect amplitude visibly the same object seen twice: in WY2015 the lines are closely spaced (shallow
lapse rate, 3.8 d/100 m) yet bowed (amplitude 7.8 d).

**Use.** This is the per-year diagnostic behind sections 12–13 and the natural supplementary sweep (one grid per
range); for the main text the (range, year) figure is the way to show a named year, e.g. WY2015 in the Sierra Nevada
next to a "usual" year, in place of describing it.
