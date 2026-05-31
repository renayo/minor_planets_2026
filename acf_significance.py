"""
acf_significance.py
    — Significance tests for autocorrelation functions (ACFs) at the
      astrological lags, using corrected longitude geometry.

  For each special point's 360-point mean-article-count wave, compute
  the biased *circular* ACF at lags corresponding to the astrological
  aspects (semi-sextile 30°, decile 36°, novile 40°, semi-square 45°,
  septile 51°, sextile 60°, quintile 72°, square 90°, trine 120°,
  opposition 180°, plus 24° quindecile). The circular estimator is the
  correct form for the periodic 360° angular wave (degree 359 wraps to
  degree 0); it matches the estimator used for the per-SP ACF figures.

  Reference 95% CI threshold for white-noise null: ±1.96/√360 ≈ 0.1033.

  Per-SP and OMNIBUS tests via:
    • K_acf       : # of astrological lags where |ACF(lag)| > 0.1033
    • joint_acf_sq: Σ ACF² over astrological lags
    • LB_Q        : Ljung-Box Q restricted to astrological lags
                    Q = n(n+2) · Σ ACF(k)² / (n − k)
                    Under iid white-noise H0, Q ~ χ²(|lags|)

  Two MC nulls:
    1. Compound (cross-planet swap + per-planet phase shift)
    2. Unnamed-MBA comparison group
  Both use corrected SP geometry and 1 122 Horizons-verified planets.

  Per-SP results report frequentist p (Bonferroni-corrected) + BF.
  Omnibus tests combine across all 11 SPs via:
    • SUM / MAX aggregates (MC-based)
    • Fisher and Stouffer combinations
"""

import csv as csv_mod
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.special import ndtri, erf
from scipy.stats import chi2 as chi2_dist
from scipy.stats import beta as beta_dist

warnings.filterwarnings('ignore')

DATA_DIR  = Path(__file__).parent
OUT_DIR   = Path(__file__).parent
CSV_FILE  = DATA_DIR / "MinorPlanetSunData.csv"
SP_NPY    = OUT_DIR / "corrected_sp_series.npy"
NAMED_NPY = OUT_DIR / "named_planet_angles_full.npy"
UNNAMED_NPY = DATA_DIR / "unnamed_mba_full_angles.npy"

N_WAVE    = 360                       # 360 integer-degree bins
N_ITER    = 5_000
SEED      = 42

# Astrological lags in degrees (= ACF lag for a 360-bin wave)
# Note: lag 0 (conjunction) is trivially 1.0 and excluded.
ASTRO_LAGS = [24, 30, 36, 40, 45, 51, 60, 72, 90, 120, 180]
# Periods: quindecile, semi-sextile, decile, novile, semi-square,
#          septile (~51), sextile, quintile, square, trine, opposition

ACF_THRESHOLD = 1.96 / np.sqrt(N_WAVE)   # ≈ 0.1033, the 95% CI for n=360

LABELS = [
    "Mesarthim", "Sun", "Moon", "Mercury", "Venus",
    "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Rahu",
]


def _parse(s):
    try:    return float(s.strip())
    except ValueError: return float("nan")


def load_csv(path):
    with open(path, newline="") as f:
        rows = list(csv_mod.reader(f))[:-1]
    return np.array([[_parse(c) for c in r[1:]] for r in rows[1:]], dtype=float)


def normalize_l2(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v.copy()


def wave_from_pairs(angle_flat, value_flat, n_deg=N_WAVE):
    sums = np.bincount(angle_flat, weights=value_flat, minlength=n_deg)
    cnts = np.bincount(angle_flat,                     minlength=n_deg).clip(1)
    return sums / cnts


def biased_acf(wave, lags):
    """Biased *circular* sample ACF at the given lags.
        ACF(k) = Σ_t (x[t] - μ)(x[(t+k) mod N] - μ) / Σ_t (x[t] - μ)²
       The full circular ACF is obtained in one step via the
       Wiener-Khinchin relation (ACF = IFFT(|FFT(x)|²) / var), then
       indexed at the requested lags. The circular form is correct for
       the periodic 360-bin angular wave, in which degree 359 is
       adjacent to degree 0; the earlier linear estimator silently
       omitted the wrap-around products and diverged from the circular
       value increasingly at the larger lags (e.g. it returned exactly
       half the circular value at the 180° opposition).
       Returns array of ACF(k) for k in lags."""
    x = np.asarray(wave, dtype=np.float64) - wave.mean()
    var = (x * x).sum()
    if var == 0.0:
        return np.zeros(len(lags))
    F = np.fft.rfft(x)
    acf_full = np.fft.irfft(F * np.conj(F), n=len(x)) / var
    return acf_full[np.asarray(lags, dtype=int)]


def stat_K_acf(acf_vals, threshold=ACF_THRESHOLD):
    """# of lags where |ACF| > threshold."""
    return int((np.abs(acf_vals) > threshold).sum())


def stat_joint_acf_sq(acf_vals):
    """Sum of squared ACF values at the astrological lags."""
    return float((acf_vals ** 2).sum())


def stat_LB_Q(acf_vals, n_obs, lags):
    """Ljung-Box Q restricted to given lags."""
    Q = n_obs * (n_obs + 2) * np.sum(acf_vals ** 2 / (n_obs - np.asarray(lags)))
    return float(Q)


def clopper_pearson_ci(rank, N, conf=0.95):
    alpha = 1 - conf
    lo = beta_dist.ppf(alpha / 2, rank, N - rank + 1) if rank > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2, rank + 1, N - rank) if rank < N else 1.0
    return float(lo), float(hi)


def interpret_bf(bf):
    if bf >= 100:    return "decisive (M_coupled)"
    if bf >= 30:     return "very strong (M_coupled)"
    if bf >= 10:     return "strong (M_coupled)"
    if bf >= 3:      return "moderate (M_coupled)"
    if bf >= 1:      return "anecdotal (M_coupled)"
    if bf >= 1/3:    return "anecdotal (M_random)"
    if bf >= 1/10:   return "moderate (M_random)"
    if bf >= 1/30:   return "strong (M_random)"
    if bf >= 1/100:  return "very strong (M_random)"
    return "decisive (M_random)"


def fisher_stouffer(per_sp_p):
    eps = 1e-10
    p = np.clip(per_sp_p, eps, 1 - eps)
    X2 = -2 * np.log(p).sum()
    df = 2 * len(p)
    p_fisher = float(1 - chi2_dist.cdf(X2, df))
    z = ndtri(1 - p)
    Z = z.sum() / np.sqrt(len(p))
    p_stouffer = float(1 - 0.5 * (1 + erf(Z / np.sqrt(2))))
    return X2, p_fisher, Z, p_stouffer


def main():
    rng = np.random.default_rng(SEED)
    t0  = time.time()

    print("=" * 70)
    print("ACF SIGNIFICANCE TESTS  —  per-SP and omnibus, corrected geometry")
    print("=" * 70)
    print(f"\nAstrological ACF lags: {ASTRO_LAGS}")
    print(f"95%-CI threshold (n=360 white noise): ±{ACF_THRESHOLD:.4f}")

    body = load_csv(CSV_FILE)
    articles_raw = body[:, 0::2].T
    n_planets_full, n_dates = articles_raw.shape

    named_angles_full = np.load(NAMED_NPY)
    valid_mask = ~np.any(np.isnan(named_angles_full), axis=1)
    n_planets = int(valid_mask.sum())
    print(f"Horizons-verified planets: {n_planets}/{n_planets_full}")

    named_angles = named_angles_full[valid_mask]
    articles = articles_raw[valid_mask]
    sp_series_corrected = np.load(SP_NPY)
    n_sp = sp_series_corrected.shape[0]

    norm_articles = np.array(
        [normalize_l2(articles[i].astype(float)) for i in range(n_planets)],
        dtype=np.float32
    )

    # Build angle tables
    print(f"\nBuilding angle tables …")
    just_angles_named = np.empty((n_sp, n_planets, n_dates), dtype=np.int16)
    for j in range(n_sp):
        sp = sp_series_corrected[j]
        ang = (named_angles - sp[None, :]) % 360.0
        just_angles_named[j] = (np.ceil(ang).astype(np.int32) % 360).astype(np.int16)
    ang_flats = [just_angles_named[j].reshape(-1).astype(np.int64) for j in range(n_sp)]

    # ── Observed ACFs ─────────────────────────────────────────────── #
    print(f"\nObserved ACF at astrological lags per SP:")
    obs_acf       = np.zeros((n_sp, len(ASTRO_LAGS)))
    obs_K_acf     = np.zeros(n_sp, dtype=int)
    obs_joint_acf = np.zeros(n_sp)
    obs_LB        = np.zeros(n_sp)
    for j in range(n_sp):
        wave = wave_from_pairs(ang_flats[j], norm_articles.reshape(-1))
        acf  = biased_acf(wave, ASTRO_LAGS)
        obs_acf[j] = acf
        obs_K_acf[j]     = stat_K_acf(acf)
        obs_joint_acf[j] = stat_joint_acf_sq(acf)
        obs_LB[j]        = stat_LB_Q(acf, N_WAVE, ASTRO_LAGS)

    # Display ACF table
    print(f"\n  {'SP':<14}" + "".join(f"  lag {l:>3d}" for l in ASTRO_LAGS) +
          f"   K_acf  joint_acf_sq  LB_Q")
    print("  " + "-" * (14 + 12 * len(ASTRO_LAGS) + 30))
    for j in range(n_sp):
        row = f"  {LABELS[j]:<14}"
        for k_idx in range(len(ASTRO_LAGS)):
            val = obs_acf[j, k_idx]
            mark = "*" if abs(val) > ACF_THRESHOLD else " "
            row += f"  {val:+6.3f}{mark}"
        row += f"   {obs_K_acf[j]:>5d}      {obs_joint_acf[j]:.4f}    {obs_LB[j]:.2f}"
        print(row)
    print(f"  (*) marks lags where |ACF| > {ACF_THRESHOLD:.4f} (95% CI for white-noise H0)")

    # ── MC nulls ──────────────────────────────────────────────────── #
    def run_mc(null_kind, just_angles_unnamed=None):
        null_K       = np.zeros((N_ITER, n_sp), dtype=np.int8)
        null_joint   = np.zeros((N_ITER, n_sp), dtype=np.float64)
        null_LB      = np.zeros((N_ITER, n_sp), dtype=np.float64)
        row_idx   = np.arange(n_planets)[:, None]
        base_cols = np.arange(n_dates)[None, :]
        rep = max(1, N_ITER // 10)
        t_start = time.time()
        for i in range(N_ITER):
            if i % rep == 0:
                elapsed = time.time() - t_start
                eta     = elapsed / max(i, 1) * (N_ITER - i)
                print(f"    iter {i:6,d}/{N_ITER:,}  elapsed {elapsed:4.0f}s  ETA {eta:4.0f}s")
            if null_kind == "compound":
                perm = rng.permutation(n_planets)
                swapped = norm_articles[perm]
                shifts = rng.integers(0, n_dates, size=n_planets)
                col_idx = (base_cols + shifts[:, None]) % n_dates
                shifted_flat = swapped[row_idx, col_idx].reshape(-1)
                for j in range(n_sp):
                    wave = wave_from_pairs(ang_flats[j], shifted_flat)
                    acf = biased_acf(wave, ASTRO_LAGS)
                    null_K[i, j]     = stat_K_acf(acf)
                    null_joint[i, j] = stat_joint_acf_sq(acf)
                    null_LB[i, j]    = stat_LB_Q(acf, N_WAVE, ASTRO_LAGS)
            else:
                perm = rng.permutation(n_planets)
                for j in range(n_sp):
                    ang_flat = just_angles_unnamed[j, perm].reshape(-1).astype(np.int64)
                    wave = wave_from_pairs(ang_flat, norm_articles.reshape(-1))
                    acf = biased_acf(wave, ASTRO_LAGS)
                    null_K[i, j]     = stat_K_acf(acf)
                    null_joint[i, j] = stat_joint_acf_sq(acf)
                    null_LB[i, j]    = stat_LB_Q(acf, N_WAVE, ASTRO_LAGS)
        print(f"  {null_kind} MC complete in {time.time() - t_start:.1f}s")
        return null_K, null_joint, null_LB

    print(f"\n[1/2] Compound-null MC …")
    null_K_cmp, null_joint_cmp, null_LB_cmp = run_mc("compound")

    print(f"\n[2/2] Unnamed-MBA-null MC …")
    unnamed_angles_full = np.load(UNNAMED_NPY)
    n_use = min(unnamed_angles_full.shape[0], n_planets)
    unnamed_angles = unnamed_angles_full[:n_use]
    just_angles_unnamed = np.empty((n_sp, n_use, n_dates), dtype=np.int16)
    for j in range(n_sp):
        sp = sp_series_corrected[j]
        ang = (unnamed_angles - sp[None, :]) % 360.0
        just_angles_unnamed[j] = (np.ceil(ang).astype(np.int32) % 360).astype(np.int16)
    null_K_mba, null_joint_mba, null_LB_mba = run_mc("unnamed_mba", just_angles_unnamed)

    # ── Per-SP stats ──────────────────────────────────────────────── #
    def per_sp(null_arr, obs_arr):
        rank      = (null_arr < obs_arr[None, :]).sum(axis=0)
        n_ge_obs  = (null_arr >= obs_arr[None, :]).sum(axis=0)
        null_mean = null_arr.mean(axis=0)
        pct       = rank / N_ITER
        p_value   = n_ge_obs / N_ITER
        p_bonf    = np.minimum(p_value * n_sp, 1.0)
        eps     = 1.0 / (N_ITER + 1)
        p_clip  = np.clip(p_value, eps, 1 - eps)
        BF      = 0.5 / p_clip
        return dict(rank=rank, p_value=p_value, p_bonf=p_bonf, BF=BF,
                    null_mean=null_mean, pct=pct)

    def stars(pb):
        if pb < 0.001: return "*** p<.001 Bonf."
        if pb < 0.01:  return "**  p<.01  Bonf."
        if pb < 0.05:  return "*   p<.05  Bonf."
        return "    n.s."

    print(f"\n\n{'='*120}")
    print("PER-SP RESULTS (ACF)")
    print('='*120)
    for null_label, nullK, nullJ, nullLB in [
        ("COMPOUND NULL", null_K_cmp, null_joint_cmp, null_LB_cmp),
        ("UNNAMED-MBA NULL", null_K_mba, null_joint_mba, null_LB_mba),
    ]:
        print(f"\n  ── {null_label} ──")
        for stat_lbl, null_arr, obs_arr in [
            ("K_acf (# astro lags >95% CI)", nullK, obs_K_acf.astype(float)),
            ("joint_acf_sq (Σ ACF²)",         nullJ, obs_joint_acf),
            ("LB_Q (Ljung-Box at astro lags)", nullLB, obs_LB),
        ]:
            r = per_sp(null_arr, obs_arr)
            print(f"\n    [{stat_lbl}]")
            print(f"    {'SP':<14} {'T_obs':>10} {'null mean':>10} {'p':>8} {'p_Bonf':>8} {'BF':>10}  Decision")
            print("    " + "-" * 90)
            for j in range(n_sp):
                print(f"    {LABELS[j]:<14} {obs_arr[j]:>10.4f} {r['null_mean'][j]:>10.4f} "
                      f"{r['p_value'][j]:>8.4f} {r['p_bonf'][j]:>8.4f} {r['BF'][j]:>10.3g}  "
                      f"{stars(r['p_bonf'][j])} ({interpret_bf(r['BF'][j])})")

    # ── OMNIBUS ──────────────────────────────────────────────────── #
    print(f"\n\n{'='*120}")
    print("OMNIBUS TESTS  (combine across all 11 SPs)")
    print('='*120)

    for null_label, nullK, nullJ, nullLB in [
        ("COMPOUND NULL", null_K_cmp, null_joint_cmp, null_LB_cmp),
        ("UNNAMED-MBA NULL", null_K_mba, null_joint_mba, null_LB_mba),
    ]:
        print(f"\n  ── {null_label} ──")
        for stat_lbl, null_arr, obs_arr in [
            ("K_acf",        nullK, obs_K_acf.astype(float)),
            ("joint_acf_sq", nullJ, obs_joint_acf),
            ("LB_Q",         nullLB, obs_LB),
        ]:
            # MC-based aggregates
            obs_sum = float(obs_arr.sum()); obs_max = float(obs_arr.max())
            null_sum = null_arr.sum(axis=1); null_max = null_arr.max(axis=1)
            p_sum = float((null_sum >= obs_sum).mean())
            p_max = float((null_max >= obs_max).mean())
            bf_sum = 0.5 / max(p_sum, 1e-5)
            bf_max = 0.5 / max(p_max, 1e-5)
            # Per-SP p for Fisher/Stouffer
            per_sp_p = (null_arr >= obs_arr[None, :]).mean(axis=0)
            X2, p_fish, Z, p_stouf = fisher_stouffer(per_sp_p)
            print(f"\n    [{stat_lbl}]")
            print(f"      MC SUM-T: obs={obs_sum:.4f}  null mean={null_sum.mean():.4f}  "
                  f"p={p_sum:.4f}  BF={bf_sum:.3g}  ({interpret_bf(bf_sum)})")
            print(f"      MC MAX-T: obs={obs_max:.4f}  null mean={null_max.mean():.4f}  "
                  f"p={p_max:.4f}  BF={bf_max:.3g}  ({interpret_bf(bf_max)})")
            print(f"      Fisher:   X²={X2:>6.2f}  p={p_fish:.6f}  "
                  f"({'sig' if p_fish < 0.05 else 'n.s.'} at α=0.05)")
            print(f"      Stouffer: Z={Z:>+6.3f}  p={p_stouf:.6f}  "
                  f"({'sig' if p_stouf < 0.05 else 'n.s.'} at α=0.05)")

    # Save CSVs
    out_csv = OUT_DIR / "acf_significance_results.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv_mod.writer(f)
        w.writerow(["special_point", "K_acf_obs", "joint_acf_sq_obs", "LB_Q_obs"] +
                   [f"ACF_lag_{l}" for l in ASTRO_LAGS])
        for j in range(n_sp):
            row = [LABELS[j], int(obs_K_acf[j]), f"{obs_joint_acf[j]:.6f}",
                   f"{obs_LB[j]:.4f}"]
            for k_idx in range(len(ASTRO_LAGS)):
                row.append(f"{obs_acf[j, k_idx]:+.6f}")
            w.writerow(row)
    print(f"\n  Saved CSV → {out_csv}")
    print(f"\nTotal runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
