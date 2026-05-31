## Effect Size Analysis

> Companion to `ACF Analysis.md`. Where `acf_significance.py` answers *is there
> an effect, and how strong is the evidence* (p-values, Bayes factors),
> `acf_effect_size.py` answers the separate question *how big is the effect*.
> The numbers below are from the 2,000-replicate run of `acf_effect_size.py`
> against both null models; `acf_effect_size_results.csv` and
> `acf_effect_size_per_lag.csv` accompany it.

**The autocorrelation effect is small in absolute correlation terms but large
relative to chance, and roughly three-quarters of it survives the strict
unnamed-asteroid control.**

### Two scales of effect size

The analysis reports the effect on two scales, which must not be conflated:

- **RMS autocorrelation** — the root-mean-square ACF across the eleven
  astrological lags, √(mean ACF²). Because an autocorrelation is itself a
  correlation coefficient, this is a bounded, absolute effect size on the
  familiar *r* scale. It answers "how strong is the typical aspect-lag
  correlation."
- **Standardized effect d** — (T_obs − mean_null) / SD_null for the
  `joint_acf_sq` statistic. It expresses how far the observed value sits from
  the null distribution, in null standard deviations, and answers "how far
  from chance." It does not grow with the number of Monte-Carlo iterations, so
  it is a genuine effect size rather than a test statistic.

A small absolute effect can sit many standard deviations from a tight null.
That is precisely the pattern here.

### Observed effect size — RMS autocorrelation

| SP          | joint_acf_sq | RMS-ACF | Magnitude  |
| ----------- | ------------ | ------- | ---------- |
| Mesarthim   | 0.546        | 0.223   | small      |
| Pluto       | 0.499        | 0.213   | small      |
| Neptune     | 0.486        | 0.210   | small      |
| Uranus      | 0.390        | 0.188   | small      |
| Saturn      | 0.383        | 0.187   | small      |
| Jupiter     | 0.317        | 0.170   | small      |
| Sun         | 0.153        | 0.118   | small      |
| Mars        | 0.136        | 0.111   | small      |
| Mercury     | 0.095        | 0.093   | negligible |
| Venus       | 0.058        | 0.073   | negligible |
| Moon        | 0.026        | 0.049   | negligible |
| **All (121)** | **3.090**  | **0.160** | **small** |

The overall RMS autocorrelation is 0.16 — "small" by Cohen's rough conventions
(≈0.1 small, 0.3 medium, 0.5 large). For reference, pure chance in this design
produces an RMS of about 0.089 under the compound null and 0.096 under the
unnamed-asteroid null, so the observed 0.16 is roughly 1.7–1.8 times the chance
level. Mesarthim, Pluto and Neptune are the strongest at ≈0.21–0.22; the Moon,
at 0.05, falls below even the null floor.

One caveat on the RMS: it averages over all eleven lags, including the
near-zero square (90°) and the weaker aspects, which pulls the figure down. At
the principal aspect lags the individual ACF values are considerably larger —
for the strong reference points the ACF reaches roughly +0.27 to +0.33 at the
small aspects (24°–36°) and about −0.25 to −0.37 at the opposition (180°), i.e.
"medium" to "large" on the *r* scale. The honest one-line characterisation is
therefore: medium-sized at the principal aspect lags, diluted to small in the
eleven-lag RMS. Per-lag values with Fisher-z 95% confidence intervals (roughly
±0.10) are in `acf_effect_size_per_lag.csv`.

### Standardized effect — distance from the null

| SP          | d (compound) | ratio | d (unnamed-MBA) | ratio |
| ----------- | ------------ | ----- | --------------- | ----- |
| Mesarthim   | 4.87         | 4.94  | 3.49            | 4.19  |
| Pluto       | 4.41         | 4.61  | 3.14            | 3.86  |
| Neptune     | 4.16         | 4.34  | 3.01            | 3.76  |
| Saturn      | 3.03         | 3.50  | 2.33            | 3.09  |
| Uranus      | 2.82         | 3.35  | 2.15            | 2.94  |
| Jupiter     | 2.23         | 2.84  | 1.81            | 2.56  |
| Sun         | 1.48         | 2.43  | 1.17            | 2.20  |
| Mercury     | 1.42         | 1.94  | 1.16            | 1.78  |
| Venus       | 0.69         | 1.42  | 0.11            | 1.09  |
| Mars        | 0.39         | 1.42  | 0.10            | 1.12  |
| Moon        | −0.59        | 0.72  | −0.64           | 0.71  |
| **Omnibus sum** | **3.28** | **3.24** | **2.39**     | **2.80** |

Measured against the null distributions, the picture is far stronger than
"small" suggests, because the null distribution of `joint_acf_sq` is narrow.
Under the compound null the strong reference points sit 3 to 5 null standard
deviations above chance (Mesarthim 4.9, Pluto 4.4, Neptune 4.2, Saturn 3.0),
and the omnibus aggregate sits 3.3 SDs above chance, with the observed
`joint_acf_sq` 3.2 times its null expectation. The Moon's standardized effect
is negative under both nulls — its observed value lies *below* the null mean —
which makes it a useful internal negative control: a reference point expected
to carry no aspect structure, and it carries none.

The omnibus *d* of 3.3 is lower than the largest per-point values because the
eleven reference points are positively correlated; summing them inflates the
null SD as well as the null mean. It is the standardized effect for the
aggregate statistic, not an average of the per-point effects.

### The named-versus-unnamed contrast

The standardized effect against the unnamed-MBA null is the contrast of central
interest: it measures how far the *named* catalogue sits above a baseline of
unnamed main-belt asteroids carrying the same kind of orbital geometry.
Comparing the two *d* columns isolates the part of the effect the names
contribute.

The effect attenuates but does not collapse. The omnibus standardized effect
falls from 3.28 to 2.39 — a reduction of about 27% — so roughly three-quarters
of the catalogue-level effect is specific to bodies that bear names, and about
one quarter is attributable to orbital geometry common to any minor planet. The
strong reference points behave the same way individually: Mesarthim, Pluto,
Neptune, Saturn, Uranus, Jupiter, the Sun and Mercury each retain 70–81% of
their standardized effect under the strict null, all still sitting 1.2 to 3.5
SDs above the unnamed baseline.

Two reference points behave differently. Mars and Venus, already weak under the
compound null (*d* ≈ 0.4 and 0.7), fall essentially to the unnamed baseline
(*d* ≈ 0.10 and 0.11) — for these two, what little departure existed is
explained almost entirely by shared geometry, with no detectable name-specific
residual. The Moon stays below the null under both models. The name-specific
signal is thus carried by the slow outer planets, Saturn, Jupiter, the Sun,
Mercury and the fixed star Mesarthim, and not by the fast inner bodies.

### Reconciling the two scales

There is no contradiction between "small" and "3.3 SDs from chance." The
absolute autocorrelation is genuinely modest — a typical |ACF| of 0.16,
individual reference points peaking near 0.22 on the RMS scale, none reaching
the medium threshold of 0.30 on that measure. But the modest pattern is
*consistent*: the same positive-at-small-aspects, negative-at-large-aspects
shape repeats across eleven lags and eight of eleven reference points, and a
small but coherent pattern repeated that many times lands far out in a tight
null distribution. The effect-size analysis and the significance analysis
therefore agree — the structure is real and robustly above chance while
remaining small in absolute magnitude. That is the appropriate expectation for
a subtle signal in media language, and a reason to report effect sizes
alongside the Bayes factors rather than in place of them.

### For the paper

> "Effect sizes were computed to complement the significance tests. On an
> absolute scale, the root-mean-square autocorrelation across the eleven
> astrological lags was 0.16 overall (Cohen's *r*: small), rising to 0.21–0.22
> for the strongest reference points and approximately 1.7–1.8 times the
> Monte-Carlo chance level; individual ACF values at the principal aspect lags
> reached the medium-to-large range. On a standardized scale, the observed
> joint autocorrelation statistic lay 3.3 null standard deviations above the
> compound null and 2.4 above the stricter unnamed-asteroid null — an
> attenuation of roughly 27% — corresponding to 3.2 and 2.8 times the
> respective null expectations. Approximately three-quarters of the
> catalogue-level effect was therefore specific to named bodies. The signal was
> carried by the slow-moving outer planets, the Sun and the fixed star
> Mesarthim; the Moon, with a negative standardized effect under both nulls,
> served as an internal negative control."

#### See [main repository](https://github.com/renayo/cosmic_semiotics/tree/main) for the analysis code and figures.
