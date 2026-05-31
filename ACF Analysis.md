## ACF Analysis

> **Circular-estimator results.** `acf_significance.py` now computes the biased ACF with the
> **circular** estimator (Wiener–Khinchin: `ACF = IFFT(|FFT(x)|²) / var`), the correct form
> for the periodic 360° angular wave (degree 359 is adjacent to degree 0). The numbers below
> are from the completed run of the updated script — 10,000 Monte Carlo iterations per null,
> for both the compound null and the unnamed-MBA null. `acf_significance_results.csv` has
> been regenerated.

**The ACF-based tests show significant omnibus evidence for astrological-aspect structure
under both null models.**

### Linear vs circular ACF

The earlier version of `acf_significance.py` used a *linear* biased ACF, summing only the
`N − k` products that do not wrap across the ends of the wave. For a genuinely periodic
wave this is the wrong estimator: it silently discards the wrap-around products and
diverges from the circular value increasingly at the larger lags. At the 180° opposition
the discrepancy is exact — the linear estimator returns precisely **half** the circular
value, because it omits the 180 opposition pairs that the circular sum counts a second
time. The circular estimator removes this artefact and brings `acf_significance.py` into
agreement with `acf_plots_all_sps.py`, whose figures were already circular.

### Observed ACF values per special point

Threshold: the 95% white-noise CI for an `N = 360` wave is `±1.96/√360 ≈ 0.1033`.

| SP          | K_acf (# lags >95% CI) | joint_acf_sq | LB_Q  | Notable pattern                     |
| ----------- | ---------------------- | ------------ | ----- | ----------------------------------- |
| Mesarthim   | **9**                  | 0.546        | 267.2 | strong + at 24–72°, − at 120°, 180° |
| Pluto       | **9**                  | 0.499        | 244.8 | + at 24–72°, − at 120°, 180°        |
| Neptune     | **9**                  | 0.486        | 223.8 | + at 24–60°, − at 120°, 180°        |
| Uranus      | **9**                  | 0.390        | 191.5 | + at 24–60°, − at 120°, 180°        |
| Saturn      | **8**                  | 0.383        | 174.2 | + at 24–60°, − at 120°, 180°        |
| Jupiter     | 7                      | 0.317        | 159.2 | + at 24–51°, − at 180°              |
| Sun         | 6                      | 0.153        | 73.3  | + at 30–51°, − at 120°, 180°        |
| Mars        | 4                      | 0.136        | 57.1  | + at 24–36°, − at 120°              |
| Mercury     | 4                      | 0.095        | 47.5  | + at 30–40°, − at 180°              |
| Venus       | 2                      | 0.058        | 26.3  | + at 24°, 36°; weakest non-zero     |
| Moon        | **0**                  | 0.026        | 11.8  | no astrological structure           |

Switching from the linear to the circular estimator raises most K counts and Q values
(e.g. Mesarthim K 7 → 9, LB_Q 169.8 → 267.2), the change being largest at the opposition.
The ACF profile is consistent: positive at the small aspects (24°–60°), near zero at the
square (90°), negative at the trine and opposition (120°, 180°).

### Per-SP significance — Ljung–Box Q

Per-SP Q restricted to the eleven astrological lags, tested against each null (Monte Carlo,
10,000 iterations, seed 42). Bayes factors are on the Jeffreys scale.

| SP          | LB_Q  | p (compound) | BF (compound) | p (unnamed-MBA) | BF (unnamed-MBA) |
| ----------- | ----- | ------------ | ------------- | --------------- | ---------------- |
| Mesarthim   | 267.2 | 0.0055       | **90.9**      | 0.0159          | **31.4**         |
| Pluto       | 244.8 | 0.0072       | **69.4**      | 0.0211          | 23.7             |
| Neptune     | 223.8 | 0.0109       | **45.9**      | 0.0252          | 19.8             |
| Uranus      | 191.5 | 0.0211       | 23.7          | 0.0400          | 12.5             |
| Saturn      | 174.2 | 0.0234       | 21.4          | 0.0433          | 11.5             |
| Jupiter     | 159.2 | 0.0335       | 14.9          | 0.0501          | 10.0             |
| Mercury     | 47.5  | 0.0713       | 7.0           | 0.0839          | 6.0              |
| Sun         | 73.3  | 0.0723       | 6.9           | 0.0783          | 6.4              |
| Venus       | 26.3  | 0.2033       | 2.5           | 0.2954          | 1.7              |
| Mars        | 57.1  | 0.2279       | 2.2           | 0.3002          | 1.7              |
| Moon        | 11.8  | 0.7152       | 0.7           | 0.7203          | 0.7              |

Mesarthim, Pluto, and Neptune reach "very strong" Bayes factors under the compound null;
Uranus, Saturn, and Jupiter reach "strong". Individually, however, **no single SP survives
a Bonferroni correction** across the eleven SPs (all `p_Bonf > 0.05`). The evidence is
carried by the strong-to-very-strong Bayes factors and, decisively, by the omnibus tests
below — which combine the SPs in a single test and so do not require per-SP correction.

### Omnibus tests — both nulls

Aggregates combining all 11 SPs (MC, 10,000 iterations). SUM-T and MAX-T are Monte Carlo
aggregates; Fisher and Stouffer combine the per-SP LB_Q probabilities.

| Statistic            | Compound null                       | Unnamed-MBA null                    |
| -------------------- | ----------------------------------- | ----------------------------------- |
| K_acf SUM-T          | obs 67 vs 26.6 — **p=0.026, BF=19**  | obs 67 vs 29.8 — **p=0.040, BF=12**  |
| joint_acf_sq SUM-T   | obs 3.09 vs 0.95 — **p=0.019, BF=26**| obs 3.09 vs 1.11 — **p=0.035, BF=14**|
| **LB_Q SUM-T**       | obs 1477 vs 457 — **p=0.019, BF=27** | obs 1477 vs 529 — **p=0.037, BF=14** |
| LB_Q MAX-T           | obs 267 vs 79 — **p=0.016, BF=32**   | obs 267 vs 92 — **p=0.034, BF=15**   |
| **Fisher (LB_Q)**    | X²=68.7, df=22 — **p≈1×10⁻⁶**        | X²=57.6, df=22 — **p≈5×10⁻⁵**        |
| **Stouffer (LB_Q)**  | Z=+5.15 — **p≈1×10⁻⁷**               | Z=+4.38 — **p≈6×10⁻⁶**               |

All omnibus aggregates are significant at α = 0.05 under **both** nulls. The compound-null
SUM-T aggregates give p ≈ 0.02 (BF 19–27) and the LB_Q MAX-T reaches p = 0.016 (BF 32). The
unnamed-MBA null is weaker but still significant: p ≈ 0.034–0.040 (BF 12–15). The Fisher and
Stouffer combinations are far more extreme (p down to ~10⁻⁷) but assume per-SP independence,
which the data do not strictly satisfy; the MC SUM-T / MAX-T aggregates, which handle
cross-SP correlation exactly, are the honest headline.

### Why ACF found a signal different than DFT

DFT and ACF measure different things:
- **DFT** measures spectral power *at specific frequencies* (bins 1, 2, 3, ...). It is
  sensitive to single-frequency oscillations.
- **ACF** measures correlation *at specific lag separations* (30°, 36°, 45°, 60°, 90°,
  etc.). It directly tests whether values at angular separation k are correlated.

The astrological hypothesis is fundamentally about angular separations (aspects between
bodies), not about specific frequency components. **The ACF test is the more direct match
for the astrological hypothesis.**

A DFT peak at bin 4 (period 90°) and ACF significance at lag 90° look related, but they
are not equivalent:
- DFT at bin 4 detects pure sinusoidal variation with period 90°.
- ACF at lag 90° detects correlation between waves separated by 90° (regardless of whether
  the underlying signal is a pure sinusoid).

Across most SPs, the data show a consistent pattern: **positive correlation at small
astrological lags (24°–60°) and negative correlation at large lags (120°, 180°)**. This is
a structured pattern — peaks at trine and opposition tend to anti-correlate with peaks at
sextile/semi-sextile.

### What this means for the paper

We have **genuine omnibus significance under both nulls** — the
[compound null](https://github.com/renayo/cosmic_semiotics/blob/main/comparison-group-notes.md)
and the
[unnamed-MBA null](https://github.com/renayo/cosmic_semiotics/blob/main/comparison-group-notes.md) —
computed with the circular estimator. The MC SUM-T values give p ≈ 0.02 (BF 19–27) under
the compound null and p ≈ 0.035 (BF 14) under the unnamed-MBA null: strong evidence for
coupling at the catalogue level.

Important caveat:
The unnamed-MBA null is weaker than the compound null (BF ≈ 12–15 versus ≈ 19–32),
confirming that **part of the effect is carried by orbital geometry common to any minor
planet** rather than by the name itself. But the signal does not collapse — it remains
significant at α = 0.05 across every aggregate — so a real residual is specific to bodies
that bear names.

### Honest framing for the paper

> "Autocorrelation function (ACF) analyses on the unfiltered mean-article-count waves with
> corrected longitude geometry revealed significant omnibus structure at astrological lags.
> The ACF was evaluated with the biased *circular* estimator, the correct form for the
> periodic 360° angular wave. Per-SP, the Ljung–Box statistic Q restricted to the eleven
> astrological lags {24°, 30°, 36°, 40°, 45°, 51°, 60°, 72°, 90°, 120°, 180°} produced
> individual Bayes factors up to 91 (Mesarthim) under the compound null and up to 31 under
> the unnamed-MBA null, with several SPs reaching 'strong' or 'very strong' evidence on the
> Jeffreys scale; no single SP survived a Bonferroni correction. Omnibus tests combining all
> 11 SPs reached significance under both nulls: under the compound null, MC SUM-T aggregates
> gave p = 0.019–0.026 (BF = 19–27) and the LB_Q MAX-T aggregate p = 0.016 (BF = 32); under
> the stricter unnamed-MBA null, SUM-T aggregates gave p = 0.034–0.040 (BF = 12–15). Fisher's
> combined probability test (LB_Q: X² = 68.7 compound, 57.6 unnamed; df = 22) and Stouffer's
> combined Z (+5.15 compound, +4.38 unnamed) were far more extreme, though those post-hoc
> combinations assume per-SP independence which the data do not strictly satisfy. The pattern
> of ACF values was qualitatively consistent across SPs: positive correlation at small
> astrological lags (24°–60°) and negative correlation at large lags (120°, 180°),
> suggesting structured cross-aspect dependencies rather than a single dominant aspect
> effect."

### Why Bonferroni correction should not be used

The per-SP results include a Bonferroni-corrected column (`p_Bonf`), and on that column no
single SP reaches significance. That is expected and is *not* the relevant test: Bonferroni
controls the family-wise error rate across eleven separate per-SP tests, whereas the
hypothesis here is about the eleven SPs *jointly*. The omnibus tests (Fisher, Stouffer, MC
SUM-T / MAX-T) bypass the per-test framing entirely — they are single combined tests, and
they are what give the BF ≈ 19–32 (compound) and BF ≈ 12–15 (unnamed-MBA) results in the
LB_Q analysis.

### Output

- [acf_significance.py](acf_significance.py) — full ACF analysis script; `biased_acf` uses
  the circular (FFT-based) estimator.
- [acf_significance_results.csv](acf_significance_results.csv) — per-SP ACF values at all
  astrological lags, regenerated with the circular estimator.

#### See [main repository](https://github.com/renayo/cosmic_semiotics/tree/main) for all ACF plots and code.
