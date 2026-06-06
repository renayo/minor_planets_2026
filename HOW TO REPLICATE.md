# HOW TO REPLICATE

A step-by-step guide to reproducing the analyses behind *Oshop & Coops
(2026), Minor Planet Orbits and Astronyms in the News* using **only** the
files in [`renayo/minor_planets_2026`](https://github.com/renayo/minor_planets_2026).

---

## What this repository contains

The repository ships the complete Python reanalysis pipeline used in the
paper: the corrected-geometry circular autocorrelation, the two
Monte-Carlo null models, the effect-size analysis, the V-test
"directional concentration" landscape, and the Gauquelin comparison,
together with the input data, the intermediate angle tables, the figures
the paper embeds, and the paper itself as a `.docx`.

**Data**

- `MinorPlanetSunData.csv` (3.1 MB) — daily Google News article counts
  paired with Sun-separation for every named minor planet in the
  catalogue.
- `minor planets names to search.csv` — the 1,211-name list searched
  against Google News.
- `corrected_sp_series.npy` — the **twelve** reference points' daily
  longitudes in the Sun-relative encoding `(body − Sun) mod 360`; a
  `12 × 290` array consumed by every downstream script. Row order:
  FPOA, Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus,
  Neptune, Pluto, **Rahu** (mean north lunar node, added in this
  revision).
- `named_planet_angles_full.npy` — daily Sun-relative angular
  separations for the named minor planets; rows whose orbital geometry
  could not be Horizons-verified across the whole window are NaN-padded
  and filtered out at load time.
- `unnamed_mba_full_angles.npy` — daily Sun-relative angular separations
  for the unnamed main-belt asteroid pool used as the unnamed-MBA null
  baseline.
- `unnamed_mba_full_names.txt` — the numeric designations accompanying
  that pool, one per row.

**Scripts**

- `make_rahu_to_sun.py` — generates Rahu's Sun-relative longitude series
  from the Meeus mean-node polynomial and augments
  `corrected_sp_series.npy` to twelve rows. Provided for reproducibility;
  the `.npy` file already contains the augmented array.
- `acf_significance.py` — the circular ACF estimator and the
  Monte-Carlo significance machinery, evaluated against both nulls.
- `acf_effect_size.py` — per-point RMS autocorrelation and the
  standardised effect *d* in null standard deviations.
- `acf_plots_all_sps.py` — twelve per-reference-point figures, each a
  linear ACF beside its polar form, plus two twelve-panel overview
  grids.
- `all v test.py` — the V-test (Rayleigh against a specified direction)
  evaluated at every integer degree from 0 to 359, against both nulls.
- `plot_circular_vtest.py` — polar plots of the V-test landscape under
  both nulls.
- `all_sps_vs_gauquelin_extracted.py` — extracts the Gauquelin reference
  shape from the included GIFs and ranks each reference point's
  correlation with it.
- `effet_size_plot.py` — the lollipop chart of *d* for every reference
  point against both nulls.
- `heatmap_acf_normalized_by_joint.py` — the benefic/malefic heatmaps
  comparing autocorrelation share at each aspect lag.

**Companion writeups**

- `ACF Analysis.md` — the per-SP and omnibus autocorrelation results in
  prose.
- `Effect Size.md` — the effect-size analysis in prose.
- `comparison-group-notes.md` — the two null models explained in plain
  English.

**Paper**

- `Oshop_Coops_minor_planets_2026.docx` — the manuscript itself,
  reflecting the twelve-reference-point analysis.

**Pre-generated CSV outputs and PNG figures** are checked in so the
analyses can be inspected without re-running. The steps below explain
what each one is and which script regenerates it.

---

## Prerequisites

- **Python 3.10+** with `numpy`, `scipy`, `pandas`, `matplotlib`, and
  `Pillow` (the Gauquelin script extracts pixels from the GIFs).
- About 1 GB of free RAM for the 5,000-replicate Monte-Carlo runs.
- A single 5,000-replicate run of `acf_significance.py` takes roughly
  three to four minutes on a modern laptop; `all v test.py` runs in
  about ninety seconds.

The Monte-Carlo iteration count was set to **5,000** in this revision
(was 10,000 in the predecessor `cosmic_semiotics` repo). The qualitative
results are unchanged at this resolution; Bayes-factor magnitudes round
to slightly different values. To restore the higher-resolution run,
change `N_ITER` near the top of `acf_significance.py` and
`acf_effect_size.py` and `N_MC` near the top of `all v test.py`.

---

## Step 1. Clone the repository

```bash
git clone https://github.com/renayo/minor_planets_2026.git
cd minor_planets_2026
```

Verify the data and script files are present:

```bash
ls -1 *.csv *.npy *.py
```

If any of the following are missing, the rest of this guide will not
work: `MinorPlanetSunData.csv`, `corrected_sp_series.npy`,
`named_planet_angles_full.npy`, `unnamed_mba_full_angles.npy`,
`acf_significance.py`, `acf_effect_size.py`, `acf_plots_all_sps.py`,
`all v test.py`, `plot_circular_vtest.py`,
`all_sps_vs_gauquelin_extracted.py`, `effet_size_plot.py`,
`heatmap_acf_normalized_by_joint.py`, `make_rahu_to_sun.py`.

---

## Step 2. (Optional) Regenerate the Rahu reference series

The shipped `corrected_sp_series.npy` already includes Rahu as row 11.
If you want to regenerate it from scratch — for example to verify the
mean-node formula or to swap in the true node instead — run:

```bash
python make_rahu_to_sun.py
```

This computes Rahu's tropical longitude for each of the 290 days at
noon UTC using the Meeus mean-node polynomial (Chapter 47), differences
it from the Sun's tropical longitude (Meeus Chapter 25) to produce the
Sun-relative series, and augments `corrected_sp_series.npy` to a
`12 × 290` array. It also writes a standalone `rahu_to_sun.npy` for any
external use.

The script prints a sanity check on Rahu's total motion: the mean north
node moves retrograde at about −19.3°/yr, so the 290-day window should
show roughly −15.3° of motion. Day 0 corresponds to noon UTC on 15
February 2022; day 289 to 1 December 2022.

---

## Step 3. Tables 1, 2, and 3 — circular ACF significance

```bash
python acf_significance.py
```

This computes the biased circular autocorrelation at the eleven aspect
lags (24°, 30°, 36°, 40°, 45°, 51°, 60°, 72°, 90°, 120°, and 180°) for
each of the **twelve** reference points, then runs 5,000 Monte-Carlo
replicates under both the compound null and the unnamed-MBA null.

**Outputs:**

- `acf_significance_results.csv` — the underlying numbers for **Table 1**
  (K, joint sum of squared autocorrelations, Ljung–Box Q, and per-SP
  BF) and **Table 2** (autocorrelations at seven representative aspect
  lags).
- `acf_all_sps.csv` — every per-point autocorrelation across all 360
  lags, used internally by `acf_plots_all_sps.py`.

The omnibus statistics that **Table 3** reports — Monte-Carlo p-values
and Bayes factors for the joint Σ ACF² and Ljung–Box Q across the
twelve reference points, plus the classical Fisher and Stouffer
combinations — are printed to stdout. The paper rounds them; the raw
values are visible in the console output.

The companion narrative is `ACF Analysis.md`.

---

## Step 4. Effect size

```bash
python acf_effect_size.py
```

Computes the root-mean-square autocorrelation across the eleven aspect
lags and the standardised effect *d* against both nulls, for every
reference point including Rahu.

**Outputs:**

- `acf_effect_size_results.csv` — RMS, magnitude category, d_compound,
  d_unnamed_mba, and the ratio of observed to null mean for each null.
- `acf_effect_size_per_lag.csv` — the per-lag breakdown with Fisher-z
  confidence intervals.

The companion narrative is `Effect Size.md`.

---

## Step 5. Figure 5 (a–l) — per-reference-point autocorrelation panels

```bash
python acf_plots_all_sps.py
```

Regenerates the **twelve** per-reference-point figures, each a linear
autocorrelation panel beside its polar form, with the 95% white-noise
confidence band drawn and the aspect lags coloured by significance.

**Outputs:** `acf_fpoa.png`, `acf_sun.png`, `acf_moon.png`,
`acf_mercury.png`, `acf_venus.png`, `acf_mars.png`, `acf_jupiter.png`,
`acf_saturn.png`, `acf_uranus.png`, `acf_neptune.png`, `acf_pluto.png`,
`acf_rahu.png` — these are **Figure 5, panels (a) through (l)** in the
order they appear in the paper. Two overview grids
`acf_grid_linear.png` and `acf_grid_polar.png` are also written and are
useful as quick all-twelve summaries.

---

## Step 6. Figure 6 — effect-size lollipop

```bash
python effet_size_plot.py
```

Reads `acf_effect_size_results.csv` and renders the lollipop chart of
*d* for every reference point against both nulls, with one-sided
significance thresholds at d = 1.645, 2.326, 3.090 marked, the Moon
shown as a negative control, and an omnibus row with diamond markers.

**Output:** `effect_size_dotplot.png` — **Figure 6**, now twelve rows
including Rahu.

---

## Step 7. Figures 3 and 4 — directional V-test landscape

```bash
python "all v test.py"
python plot_circular_vtest.py
```

The first command evaluates the V-statistic at every integer degree
from 0 to 359 for each reference point, against both nulls, with 5,000
Monte-Carlo replicates per null. The second command turns the resulting
omnibus CSV into the polar plots embedded in the paper.

**Outputs:**

- `conjunction_vtest_per_sp.csv` — 4,320 rows (360 directions × 12
  reference points) with V, the standardised effect *d*, *p*, and the
  Bayes factor under both nulls.
- `conjunction_vtest_omnibus.csv` — 360 rows, one per direction, with
  the omnibus T = Σ d² combined across the twelve reference points.
- `cmp_vtest_polar_plot.png` — **Figure 3**, the compound-null
  landscape.
- `mba_vtest_polar_plot.png` — **Figure 4**, the unnamed-MBA-null
  landscape.

---

## Step 8. Figures 7 and 8 — Gauquelin comparison

```bash
python all_sps_vs_gauquelin_extracted.py
```

Extracts the classical Gauquelin reference shape from
`gauquelinzones.gif` and `gauquelinzonesrotated.gif` (both **Figure 7**),
then correlates each reference point's smoothed signal against that
shape at zero rotational shift and at the best-fitting rotation.

**Outputs:**

- `all_sps_vs_gauquelin_extracted_summary.csv` — per-SP zero-shift
  Pearson r, r², best-rotation offset, and best-rotation r.
- `all_sps_vs_gauquelin_extracted_curves.csv` — the underlying smoothed
  curves for inspection.
- `all_sps_vs_gauquelin_extracted.png` and
  `all_sps_vs_gauquelin_linear.png` — additional diagnostic views,
  twelve panels each.
- `all_sps_vs_gauquelin_rank.png` — **Figure 8**, the bar ranking of
  Pearson correlations between each reference point's residual and the
  Gauquelin reference shape, now including Rahu among the twelve.

---

## Step 9. Figure 9 and the benefic/malefic heatmaps

```bash
python heatmap_acf_normalized_by_joint.py
```

Produces the diverging-bar plot and two heatmaps that compare
autocorrelation share at each aspect lag between the classical
benefic and malefic groupings, with Rahu placed in its own
"Lunar Node (Vedic)" group at the bottom.

**Outputs:**

- `heatmap_acf_norm_by_joint_all_grouped.png` — full twelve-SP heatmap
  arranged by traditional group.
- `heatmap_acf_norm_by_joint_benefics_vs_malefics.png` — restricted
  benefic-vs-malefic comparison heatmap.
- `heatmap_acf_variance_share_benefics_vs_malefics.png` — same
  comparison expressed as variance share per row.
- `acf_norm_by_joint.csv` — the underlying normalised values.

The standalone diverging-bar variant `benefic_malefic_diverging.png` is
checked in as **Figure 9** and was produced by a separate companion
visualisation; the heatmaps above are the supersets that reach print.

---

## Step 10. Figures 1 and 2 — already in the repository

These two handbuilt figures are static and do not require regeneration:

- `wave_types.jpg` — **Figure 1**, the wave-type primer (amplitude,
  frequency, phase).
- `planetary_degrees_feb15_2022_clockwise.png` — **Figure 2**, the
  angular-difference example showing the minor planet Tyson 150
  degrees ahead of the Sun on 15 February 2022.

---

## Note on the Rahu reference point

Rahu is the mean north lunar node — a mathematical point describing the
intersection of the Moon's orbital plane with the ecliptic, central to
Vedic astrology. It moves retrograde at about −19.3°/yr, so over the
290-day study window it traverses only about fifteen degrees, joining
the slow-moving reference cluster (FPOA, Saturn, Uranus, Neptune,
Pluto) on geometric grounds.

The longitude formula used is Meeus, *Astronomical Algorithms*, 2nd
ed., Chapter 47 (mean longitude of the ascending node):

```
Ω = 125.0445550 − 1934.1361849·T + 0.0020762·T² + T³/467410 − T⁴/60616000
```

where T is Julian centuries from J2000.0. The true node — which adds a
perturbation term of up to about 1.5° — could be substituted by
swapping the implementation in `make_rahu_to_sun.py`; the Vedic
tradition uses the mean node, which is what this analysis adopts.

Rahu's results appear in Tables 1, 2, 3, and 4 alongside the other
eleven reference points; it appears as Figure 5(l) in the per-SP panel
grid, as the twelfth row in Figure 6, and as the twelfth entry in the
Figure 8 ranking.

---

## Citation

Oshop, R., & Coops, A. (2026). *Minor planet orbits and astronyms in
the news* [Data set and code]. GitHub.
<https://github.com/renayo/minor_planets_2026>
