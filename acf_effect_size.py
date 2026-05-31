"""
acf_effect_size.py
    — Effect-size measures for the autocorrelation analysis.

  acf_significance.py reports p-values and Bayes factors: it answers
  "is there an effect, and how strong is the evidence?". This companion
  script answers the separate question "how big is the effect?", with
  measures that do not grow with the sample size or with the number of
  Monte-Carlo iterations:

    - per-lag ACF          each ACF value is already a standardized effect
                           size (a correlation r); reported with a Fisher-z
                           95% confidence interval.
    - RMS autocorrelation  sqrt(mean(ACF^2)) over the 11 astrological lags
                           — one bounded effect size per special point and
                           one overall.
    - standardized effect  (T_obs - mean(T_null)) / SD(T_null) for the
                           joint_acf_sq statistic, against each null model
                           — a Monte-Carlo analogue of Cohen's d / Glass's D.
                           This does NOT shrink as iterations grow, so it is
                           a genuine effect size, not a test statistic.
    - ratio effect size    T_obs / mean(T_null), in plain units.

  The standardized effect against the unnamed-MBA null IS the
  named-vs-unnamed contrast: how far the named catalogue sits above the
  unnamed-asteroid baseline, in null-SD units.

  A note on confidence intervals: a planet-level bootstrap is deliberately
  NOT used for the RMS effect size. Resampling planets with replacement
  injects extra noise into the mean wave, and added noise dilutes a
  *normalized* autocorrelation (it inflates the variance in the denominator),
  so a bootstrap CI for ACF / RMS-ACF is biased downward. The per-lag
  Fisher-z CIs and the standardized-vs-null effect carry the uncertainty
  information instead.

  Data loading, the circular ACF estimator and the wave construction are
  imported from acf_significance.py, so every number is computed with the
  same primitives as that script. Place this file next to
  acf_significance.py and run it directly:  python acf_effect_size.py
"""

import os
import csv
import time
import numpy as np

import acf_significance as base   # reuses load_csv, normalize_l2,
                                  # wave_from_pairs, biased_acf (circular), etc.

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
N_MC = 2000        # null replicates. An effect size needs only a stable
                   # null MEAN and SD, which converge quickly; the 10,000-
                   # iteration p-value tails are acf_significance.py's job.
Z95  = 1.959964    # standard-normal 97.5th percentile (for 95% CIs)

LAGS   = np.asarray(base.ASTRO_LAGS, dtype=int)
N_LAGS = len(LAGS)
N_WAVE = base.N_WAVE
LABELS = base.LABELS


# ----------------------------------------------------------------------
# Effect-size helpers
# ----------------------------------------------------------------------
def fisher_ci(r, n, z=Z95):
    """95% CI for a correlation r via the Fisher z-transform.
       n is the wave length (360). Standard large-sample approximation."""
    r = np.clip(np.asarray(r, dtype=float), -0.999999, 0.999999)
    zt = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    return np.tanh(zt - z * se), np.tanh(zt + z * se)


def magnitude(r):
    """Cohen's rough verbal label for a correlation-type effect size."""
    a = abs(float(r))
    if a < 0.10:
        return "negligible"
    if a < 0.30:
        return "small"
    if a < 0.50:
        return "medium"
    return "large"


def acf_at_lags(wave):
    """Circular ACF at the astrological lags (imported estimator)."""
    return base.biased_acf(wave, LAGS)


# ----------------------------------------------------------------------
def main():
    rng = np.random.default_rng(base.SEED)
    print("=" * 70)
    print("ACF EFFECT SIZE  —  how big is the effect")
    print("=" * 70)

    # ---- data (loaded exactly as in acf_significance.py) ----
    body = base.load_csv(base.CSV_FILE)
    articles_raw = body[:, 0::2].T
    n_planets_full, n_dates = articles_raw.shape
    named_full = np.load(base.NAMED_NPY)
    valid = ~np.any(np.isnan(named_full), axis=1)
    named_angles = named_full[valid]
    articles = articles_raw[valid]
    n_planets = articles.shape[0]
    sp_series = np.load(base.SP_NPY)
    n_sp = sp_series.shape[0]
    norm_articles = np.array(
        [base.normalize_l2(articles[i].astype(float)) for i in range(n_planets)],
        dtype=np.float64)
    print(f"\nHorizons-verified planets : {n_planets}/{n_planets_full}")
    print(f"Astrological lags ({N_LAGS:2d})     : {LAGS.tolist()}")
    print(f"White-noise 95% CI band   : +/-{1.96/np.sqrt(N_WAVE):.4f}")

    # named-planet angle tables, shape (n_sp, n_planets, n_dates)
    ang2d = np.empty((n_sp, n_planets, n_dates), dtype=np.int64)
    for j in range(n_sp):
        ang = (named_angles - sp_series[j][None, :]) % 360.0
        ang2d[j] = np.ceil(ang).astype(np.int64) % 360
    ang_flat = [ang2d[j].reshape(-1) for j in range(n_sp)]

    # ---- observed ACF and the effect sizes that need no null model ----
    obs_acf = np.zeros((n_sp, N_LAGS))
    for j in range(n_sp):
        wave = base.wave_from_pairs(ang_flat[j], norm_articles.reshape(-1))
        obs_acf[j] = acf_at_lags(wave)
    obs_joint   = (obs_acf ** 2).sum(axis=1)            # joint_acf_sq per SP
    obs_rms     = np.sqrt((obs_acf ** 2).mean(axis=1))  # RMS-ACF per SP
    overall_rms = float(np.sqrt((obs_acf ** 2).mean())) # RMS over all 11x11
    acf_lo, acf_hi = fisher_ci(obs_acf, N_WAVE)          # Fisher-z CI per lag

    # ---- null distributions of joint_acf_sq for the standardized effects ----
    def run_null(kind, unnamed_ang=None):
        nj = np.zeros((N_MC, n_sp))
        row = np.arange(n_planets)[:, None]
        cols = np.arange(n_dates)[None, :]
        t = time.time()
        for i in range(N_MC):
            if kind == "compound":
                perm = rng.permutation(n_planets)
                sw = norm_articles[perm]
                sh = rng.integers(0, n_dates, n_planets)
                flat = sw[row, (cols + sh[:, None]) % n_dates].reshape(-1)
                for j in range(n_sp):
                    w = base.wave_from_pairs(ang_flat[j], flat)
                    nj[i, j] = (acf_at_lags(w) ** 2).sum()
            else:  # unnamed-MBA
                perm = rng.permutation(n_planets)
                flat = norm_articles.reshape(-1)
                for j in range(n_sp):
                    w = base.wave_from_pairs(
                        unnamed_ang[j][perm].reshape(-1), flat)
                    nj[i, j] = (acf_at_lags(w) ** 2).sum()
        print(f"      {kind} null done in {time.time()-t:.0f}s")
        return nj

    print(f"\n[1/2] Compound-null replicates ({N_MC}) ...")
    null_cmp = run_null("compound")

    print(f"\n[2/2] Unnamed-MBA-null replicates ({N_MC}) ...")
    if os.path.exists(base.UNNAMED_NPY):
        pool = np.load(base.UNNAMED_NPY)
        if pool.shape[0] >= n_planets:
            pool = pool[:n_planets]
            unnamed_ang = np.empty((n_sp, n_planets, n_dates), dtype=np.int64)
            for j in range(n_sp):
                ang = (pool - sp_series[j][None, :]) % 360.0
                unnamed_ang[j] = np.ceil(ang).astype(np.int64) % 360
            null_mba = run_null("unnamed_mba", unnamed_ang)
        else:
            null_mba = None
            print(f"      unnamed pool has only {pool.shape[0]} bodies "
                  f"(< {n_planets}); skipping unnamed-MBA contrast.")
    else:
        null_mba = None
        print(f"      {os.path.basename(str(base.UNNAMED_NPY))} not found "
              f"-- skipping the unnamed-MBA contrast.")

    # ---- standardized + ratio effect sizes ----
    def effects(obs, null):
        mu = null.mean(axis=0)
        sd = null.std(axis=0, ddof=1)
        return (obs - mu) / sd, obs / mu          # standardized d, ratio

    d_cmp, r_cmp = effects(obs_joint, null_cmp)
    sumO = float(obs_joint.sum())
    cmp_sum = null_cmp.sum(axis=1)
    d_cmp_sum = (sumO - cmp_sum.mean()) / cmp_sum.std(ddof=1)
    r_cmp_sum = sumO / cmp_sum.mean()

    if null_mba is not None:
        d_mba, r_mba = effects(obs_joint, null_mba)
        mba_sum = null_mba.sum(axis=1)
        d_mba_sum = (sumO - mba_sum.mean()) / mba_sum.std(ddof=1)
        r_mba_sum = sumO / mba_sum.mean()
    else:
        d_mba = np.full(n_sp, np.nan)
        r_mba = np.full(n_sp, np.nan)
        d_mba_sum = r_mba_sum = float("nan")

    def cell(x, w=8, p=3):
        return f"{x:>{w}.{p}f}" if np.isfinite(x) else f"{'n/a':>{w}}"

    # ---- report: observed effect sizes -------------------------------- #
    print("\n" + "=" * 72)
    print("OBSERVED EFFECT SIZE PER SPECIAL POINT")
    print("  RMS-ACF = typical |ACF| across the 11 aspect lags (a bounded,")
    print("  correlation-scale effect size). joint_acf_sq = sum of ACF^2.")
    print("=" * 72)
    print(f"  {'SP':<12}{'joint_acf_sq':>14}{'RMS-ACF':>10}   magnitude")
    print("  " + "-" * 68)
    order = np.argsort(-obs_rms)
    for j in order:
        print(f"  {LABELS[j]:<12}{obs_joint[j]:>14.4f}{obs_rms[j]:>10.3f}"
              f"   {magnitude(obs_rms[j])}")
    print("  " + "-" * 68)
    print(f"  {'ALL (121)':<12}{sumO:>14.4f}{overall_rms:>10.3f}"
          f"   {magnitude(overall_rms)}")
    print("\n  Per-lag ACF values (each a correlation, with Fisher-z 95% CI")
    print("  of roughly +/-0.10) are written to acf_effect_size_per_lag.csv.")

    # ---- report: standardized effects vs the nulls -------------------- #
    print("\n" + "=" * 72)
    print("STANDARDIZED EFFECT SIZE  vs the Monte-Carlo nulls  "
          "[statistic: joint_acf_sq]")
    print("  d = (T_obs - mean_null) / SD_null      ratio = T_obs / mean_null")
    print("=" * 72)
    print(f"  {'SP':<12}{'d (compound)':>14}{'ratio':>9}"
          f"{'d (unnamed)':>14}{'ratio':>9}")
    print("  " + "-" * 68)
    for j in order:
        print(f"  {LABELS[j]:<12}{cell(d_cmp[j],14)}{cell(r_cmp[j],9,2)}"
              f"{cell(d_mba[j],14)}{cell(r_mba[j],9,2)}")
    print("  " + "-" * 68)
    print(f"  {'OMNIBUS SUM':<12}{cell(d_cmp_sum,14)}{cell(r_cmp_sum,9,2)}"
          f"{cell(d_mba_sum,14)}{cell(r_mba_sum,9,2)}")

    print("\nThe standardized effect against the unnamed-MBA null is the")
    print("named-vs-unnamed contrast (how far the named catalogue sits above")
    print("the unnamed-asteroid baseline, in null-SD units). Effect size")
    print("answers 'how big'; acf_significance.py answers 'is it there'.")

    # ---- CSV: per-SP summary ------------------------------------------ #
    out1 = base.OUT_DIR / "acf_effect_size_results.csv"
    with open(out1, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["special_point", "joint_acf_sq", "rms_acf", "magnitude",
                    "d_compound", "ratio_compound",
                    "d_unnamed_mba", "ratio_unnamed_mba"])
        for j in range(n_sp):
            w.writerow([LABELS[j], f"{obs_joint[j]:.6f}", f"{obs_rms[j]:.6f}",
                        magnitude(obs_rms[j]),
                        f"{d_cmp[j]:.6f}", f"{r_cmp[j]:.6f}",
                        f"{d_mba[j]:.6f}", f"{r_mba[j]:.6f}"])
        w.writerow(["ALL_SPS", f"{sumO:.6f}", f"{overall_rms:.6f}",
                    magnitude(overall_rms),
                    f"{d_cmp_sum:.6f}", f"{r_cmp_sum:.6f}",
                    f"{d_mba_sum:.6f}", f"{r_mba_sum:.6f}"])

    # ---- CSV: per-lag ACF effect sizes with Fisher-z CIs -------------- #
    out2 = base.OUT_DIR / "acf_effect_size_per_lag.csv"
    with open(out2, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["special_point", "lag_deg", "acf",
                    "fisher_ci_lo", "fisher_ci_hi"])
        for j in range(n_sp):
            for k in range(N_LAGS):
                w.writerow([LABELS[j], int(LAGS[k]), f"{obs_acf[j, k]:.6f}",
                            f"{acf_lo[j, k]:.6f}", f"{acf_hi[j, k]:.6f}"])

    print(f"\nSaved -> {out1}")
    print(f"Saved -> {out2}")


if __name__ == "__main__":
    main()
