## ACF Analysis

> **Twelve reference points, circular-estimator results.** `acf_significance.py`
> now computes the biased ACF with the **circular** estimator (Wiener–Khinchin:
> `ACF = IFFT(|FFT(x)|²) / var`), the correct form for the periodic 360° angular
> wave (degree 359 is adjacent to degree 0). The numbers below are from the
> completed run of the updated script over **twelve** reference points — the
> original eleven (Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn, Uranus,
> Neptune, Pluto, Mesarthim) plus the mean north lunar node **Rahu**, with
> 5,000 Monte Carlo iterations per null for both the compound null and the
> unnamed-MBA null. `acf_significance_results.csv` has been regenerated.

**The ACF-based tests show significant omnibus evidence for astrological-aspect
structure under both null models, and the inclusion of Rahu strengthens both
the per-point and the omnibus signal.**

### Linear vs circular ACF

The earlier version of `acf_significance.py` used a *linear* biased ACF, summing
only the `N − k` products that do not wrap across the ends of the wave. For a
genuinely periodic wave this is the wrong estimator: it silently discards the
wrap-around products and diverges from the circular value increasingly at the
larger lags. At the 180° opposition the discrepancy is exact — the linear
estimator returns precisely **half** the circular value, because it omits the
180 opposition pairs that the circular sum counts a second time. The circular
estimator removes this artefact and brings `acf_significance.py` into
agreement with `acf_plots_all_sps.py`, whose figures were already circular.

### Observed ACF values per special point

Threshold: the 95% white-noise CI for an `N = 360` wave is `±1.96/√360 ≈ 0.1033`.

| SP        | K_acf (# lags >95% CI) | joint_acf_sq | LB_Q  | Notable pattern                     |
| --------- | ----------------------: | -----------: | ----: | ----------------------------------- |
| Mesarthim | **9**                   | 0.546        | 267.2 | strong + at 24–72°, − at 120°, 180° |
| Pluto     | **9**                   | 0.499        | 244.8 | + at 24–72°, − at 120°, 180°        |
| Neptune   | **9**                   | 0.486        | 223.8 | + at 24–60°, − at 120°, 180°        |
| Uranus    | **9**                   | 0.390        | 191.5 | + at 24–60°, − at 120°, 180°        |
| **Rahu**  | **9**                   | 0.386        | 168.6 | + at 24–60°, − at 120°, 180°        |
| Saturn    | **8**                   | 0.383        | 174.2 | + at 24–60°, − at 120°, 180°        |
| Jupiter   | 7                       | 0.317        | 159.2 | + at 24–51°, − at 180°              |
| Sun       | 6                       | 0.153        |  73.3 | + at 30–51°, − at 120°, 180°        |
| Mars      | 4                       | 0.136        |  57.1 | + at 24–36°, − at 120°              |
| Mercury   | 4                       | 0.095        |  47.5 | + at 30–40°, − at 180°              |
| Venus     | 2                       | 0.058        |  26.3 | + at 24°, 36°; weakest non-zero     |
| Moon      | **0**                   | 0.026        |  11.8 | no astrological structure           |

Rahu joins the strong-signal cluster on its own merits: nine of its eleven
astrological lags clear the white-noise threshold, the same count as Mesarthim,
Uranus, Neptune, and Pluto. Its joint sum of squared autocorrelations (0.386)
sits between Uranus (0.390) and Saturn (0.383), and its Ljung–Box Q (168.6)
sits between Saturn (174.2) and Jupiter (159.2). The qualitative ACF profile
is consistent with the slow-moving outer-body cluster: positive correlation at
the small aspects (24°–60°), near zero at the square (90°), negative at the
trine and opposition (120°, 180°). This profile is unsurprising in retrospect:
Rahu moves retrograde at about 0.05° per day, so over the 290-day window it
traverses only ~15° of ecliptic, behaving geometrically much like the slow
outer planets.

### Per-SP significance — Ljung–Box Q

Per-SP Q restricted to the eleven astrological lags, tested against each null
(Monte Carlo, 5,000 iterations, seed 42). Bayes factors are on the Jeffreys
scale.

| SP        | LB_Q  | p (compound) | BF (compound) | p (unnamed-MBA) | BF (unnamed-MBA) |
| --------- | ----: | -----------: | ------------: | --------------: | ---------------: |
| Mesarthim | 267.2 | 0.0058       | **86.2**      | 0.0186          | **26.9**         |
| Pluto     | 244.8 | 0.0068       | **73.5**      | 0.0200          | **25.0**         |
| Neptune   | 223.8 | 0.0102       | **49.0**      | 0.0272          | 18.4             |
| Uranus    | 191.5 | 0.0206       | 24.3          | 0.0402          | 12.4             |
| Saturn    | 174.2 | 0.0228       | 21.9          | 0.0446          | 11.2             |
| **Rahu**  | 168.6 | 0.0260       | **19.2**      | 0.0518          | 9.65             |
| Jupiter   | 159.2 | 0.0316       | 15.8          | 0.0538          | 9.29             |
| Sun       |  73.3 | 0.0694       | 7.20          | 0.0902          | 5.54             |
| Mars      |  57.1 | 0.2292       | 2.18          | 0.3088          | 1.62             |
| Mercury   |  47.5 | 0.0732       | 6.83          | 0.0804          | 6.22             |
| Venus     |  26.3 | 0.2040       | 2.45          | 0.3074          | 1.63             |
| Moon      |  11.8 | 0.7126       | 0.70          | 0.7282          | 0.69             |

Mesarthim, Pluto, and Neptune reach "very strong" Bayes factors under the
compound null; Uranus, Saturn, Rahu, and Jupiter reach "strong". Rahu's BF
of 19.2 sits between Saturn (21.9) and Jupiter (15.8). Under the stricter
unnamed-MBA null Rahu falls just below the Jeffreys "strong" threshold of
10 (BF = 9.65), in the same band as Jupiter (9.29). Individually, however,
**no single SP survives a Bonferroni correction** across the twelve SPs
(all `p_Bonf > 0.05`). The evidence is carried by the strong-to-very-strong
Bayes factors and, decisively, by the omnibus tests below — which combine
the SPs in a single test and so do not require per-SP correction.

### Omnibus tests — both nulls

Aggregates combining all 12 SPs (MC, 5,000 iterations). SUM-T and MAX-T are
Monte Carlo aggregates; Fisher and Stouffer combine the per-SP LB_Q
probabilities.

| Statistic              | Compound null                          | Unnamed-MBA null                       |
| ---------------------- | -------------------------------------- | -------------------------------------- |
| K_acf SUM-T            | obs 76 vs 29.6 — **p=0.024, BF=20.8**  | obs 76 vs 33.2 — **p=0.040, BF=12.4**  |
| joint_acf_sq SUM-T     | obs 3.48 vs 1.06 — **p=0.020, BF=24.8**| obs 3.48 vs 1.24 — **p=0.040, BF=12.6**|
| joint_acf_sq MAX-T     | obs 0.546 vs 0.164 — **p=0.018, BF=27.8** | obs 0.546 vs 0.192 — **p=0.036, BF=13.9** |
| **LB_Q SUM-T**         | obs 1645 vs 508 — **p=0.020, BF=24.5** | obs 1645 vs 591 — **p=0.040, BF=12.4** |
| **LB_Q MAX-T**         | obs 267 vs 80 — **p=0.016, BF=31.2**   | obs 267 vs 93 — **p=0.034, BF=14.6**   |
| **Fisher (LB_Q)**      | X²=76.4, df=24 — **p≈3×10⁻⁷**          | X²=62.6, df=24 — **p≈3×10⁻⁵**          |
| **Stouffer (LB_Q)**    | Z=+5.51 — **p≈2×10⁻⁸**                 | Z=+4.59 — **p≈2×10⁻⁶**                 |

All omnibus aggregates are significant at α = 0.05 under **both** nulls. The
compound-null SUM-T aggregates give p ≈ 0.020 (BF 20.8–24.8) and the LB_Q
MAX-T reaches p = 0.016 (BF 31.2). The unnamed-MBA null is weaker but still
significant: p ≈ 0.034–0.040 (BF 12.4–14.6). Adding Rahu to the analysis
lifted the joint_acf_sq SUM-T from 3.09 to 3.48 and the LB_Q SUM-T from
1477 to 1645, both with proportionate null-mean increases, so the p-values
and Bayes factors are comparable to the eleven-SP results — Rahu contributes
its share without disrupting the omnibus signal. The Fisher and Stouffer
combinations are far more extreme (p down to ~10⁻⁸) but assume per-SP
independence, which the data do not strictly satisfy; the MC SUM-T / MAX-T
aggregates, which handle cross-SP correlation exactly, are the honest headline.

### Why ACF found a signal different than DFT

DFT and ACF measure different things:

- **DFT** measures spectral power *at specific frequencies* (bins 1, 2, 3, …).
  It is sensitive to single-frequency oscillations.
- **ACF** measures correlation *at specific lag separations* (30°, 36°, 45°,
  60°, 90°, etc.). It directly tests whether values at angular separation k
  are correlated.

The astrological hypothesis is fundamentally about angular separations (aspects
between bodies), not about specific frequency components. **The ACF test is
the more direct match for the astrological hypothesis.**

A DFT peak at bin 4 (period 90°) and ACF significance at lag 90° look related,
but they are not equivalent:

- DFT at bin 4 detects pure sinusoidal variation with period 90°.
- ACF at lag 90° detects correlation between waves separated by 90° (regardless
  of whether the underlying signal is a pure sinusoid).

Across most SPs, the data show a consistent pattern: **positive correlation at
small astrological lags (24°–60°) and negative correlation at large lags (120°,
180°)**. This is a structured pattern — peaks at trine and opposition tend to
anti-correlate with peaks at sextile/semi-sextile.

### What this means for the paper

We have **genuine omnibus significance under both nulls** — the
[compound null](https://github.com/renayo/minor_planets_2026/blob/main/comparison-group-notes.md)
and the
[unnamed-MBA null](https://github.com/renayo/minor_planets_2026/blob/main/comparison-group-notes.md) —
computed with the circular estimator over twelve reference points. The MC
SUM-T values give p ≈ 0.020 (BF 20–28) under the compound null and p ≈ 0.040
(BF 12–14) under the unnamed-MBA null: strong evidence for coupling at the
catalogue level.

**The Rahu finding.** Adding the mean north lunar node, a Vedic-astrological
reference point with no physical body associated with it, produces a per-SP
ACF profile statistically indistinguishable from those of the slow-moving
outer planets. Nine of eleven astrological lags clear the 95% white-noise
threshold, joint Σ ACF² = 0.386, LB_Q = 168.6, BF = 19.2 (compound) and
9.65 (unnamed-MBA). The directional signal in the V-test analysis is also
nontrivial: peak BF ≈ 92 (compound) and 35 (unnamed-MBA) at direction ≈ 125°,
midway between trine and square axes. That a calculated point with no
physical body behaves like the named slow-movers is consistent with a broader
geometric interpretation of the catalogue-wide signal: what is being measured
is something about angular geometry in the tropical zodiac, not something
about the body whose name happens to be attached to the reference point.

**Important caveat.** The unnamed-MBA null is weaker than the compound null
(BF ≈ 12–15 versus ≈ 21–31), confirming that **part of the effect is carried
by orbital geometry common to any minor planet** rather than by the name
itself. But the signal does not collapse — it remains significant at α = 0.05
across every aggregate — so a real residual is specific to bodies that bear
names.

### Honest framing for the paper

> "Autocorrelation function (ACF) analyses on the unfiltered mean-article-count
> waves with corrected longitude geometry revealed significant omnibus structure
> at astrological lags. The ACF was evaluated with the biased *circular*
> estimator, the correct form for the periodic 360° angular wave. Twelve
> reference points were tested — the Sun, the Moon, the seven classical and
> modern planets, Pluto, the fixed star Mesarthim (a proxy for the First Point
> of Aries), and the mean north lunar node Rahu. Per-SP, the Ljung–Box statistic
> Q restricted to the eleven astrological lags {24°, 30°, 36°, 40°, 45°, 51°,
> 60°, 72°, 90°, 120°, 180°} produced individual Bayes factors up to 86
> (Mesarthim) under the compound null and up to 27 under the unnamed-MBA null,
> with Mesarthim, Pluto, Neptune, Uranus, Saturn, Rahu, and Jupiter all reaching
> 'strong' or 'very strong' evidence on the Jeffreys scale; no single SP
> survived a Bonferroni correction. Omnibus tests combining all 12 SPs reached
> significance under both nulls: under the compound null, MC SUM-T aggregates
> gave p = 0.020–0.024 (BF = 20.8–24.8) and the LB_Q MAX-T aggregate p = 0.016
> (BF = 31.2); under the stricter unnamed-MBA null, SUM-T aggregates gave
> p = 0.036–0.040 (BF = 12.4–13.9). Fisher's combined probability test (LB_Q:
> X² = 76.4 compound, 62.6 unnamed; df = 24) and Stouffer's combined Z (+5.51
> compound, +4.59 unnamed) were far more extreme, though those post-hoc
> combinations assume per-SP independence which the data do not strictly satisfy.
> The pattern of ACF values was qualitatively consistent across SPs: positive
> correlation at small astrological lags (24°–60°) and negative correlation at
> large lags (120°, 180°), suggesting structured cross-aspect dependencies
> rather than a single dominant aspect effect. Rahu, a mathematical point with
> no physical body, behaved geometrically like the slow-moving named bodies
> (Saturn through Pluto and Mesarthim), with per-SP BF = 19.2 (compound) and
> 9.65 (unnamed-MBA)."

### Why Bonferroni correction should not be used

The per-SP results include a Bonferroni-corrected column (`p_Bonf`), and on
that column no single SP reaches significance. That is expected and is *not*
the relevant test: Bonferroni controls the family-wise error rate across
twelve separate per-SP tests, whereas the hypothesis here is about the twelve
SPs *jointly*. The omnibus tests (Fisher, Stouffer, MC SUM-T / MAX-T) bypass
the per-test framing entirely — they are single combined tests, and they are
what give the BF ≈ 20–31 (compound) and BF ≈ 12–15 (unnamed-MBA) results in
the LB_Q analysis.

### Output

- [acf_significance.py](https://github.com/renayo/minor_planets_2026/blob/main/acf_significance.py) — full ACF analysis script; `biased_acf` uses the
  circular (FFT-based) estimator; `LABELS` now contains twelve entries
  including Rahu.
- [acf_significance_results.csv](https://github.com/renayo/minor_planets_2026/blob/main/acf_significance_results.csv) — per-SP ACF values at all
  astrological lags, regenerated with the circular estimator over twelve SPs.
- [make_rahu_to_sun.py](https://github.com/renayo/minor_planets_2026/blob/main/make_rahu_to_sun.py) — helper that generates Rahu's Sun-relative
  longitude series from the Meeus mean-node polynomial (Chapter 47) and
  augments `corrected_sp_series.npy` to twelve rows.

#### See [main repository](https://github.com/renayo/minor_planets_2026/tree/main) for all ACF plots and code.
