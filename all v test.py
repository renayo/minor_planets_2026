"""
conjunction_vtest_360.py
    — V-test for concentration across all 360 reference directions.

  This script answers it with the V-test -- the Rayleigh test of circular
  uniformity against concentration toward a *specified* direction. Here, it
  computes the V-test simultaneously for ALL 360 integer degrees (0 to 359).
  Each (body, date) pair contributes its angular separation theta from the
  reference point, weighted by that day's L2-normalised article count. The
  V-statistic is the article-weighted mean of cos(theta).

  The script generates Monte-Carlo null models and calculates the exact
  p-values, Bayes Factors, and effect sizes (d) for every single degree, 
  exporting the entire 360-degree landscape to CSV.

  Run it directly:  python conjunction_vtest_360.py
"""

import os
import csv
import time
from pathlib import Path
import numpy as np

import acf_significance as base   # wave_from_pairs, normalize_l2, load_csv,
                                  # and the constants N_WAVE / LABELS / SEED.
                                  # Path constants are NOT used.

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
N_MC      = 5000   # null replicates. A p-value test needs a well-resolved
                    # tail, so this matches acf_significance.py.

N_WAVE = base.N_WAVE          # 360
LABELS = base.LABELS

# Shape: (360 directions, 360 wave bins)
# _COS[d, a] gives the cos weight for wave bin 'a' relative to direction 'd'
_COS = np.cos(np.deg2rad(np.arange(N_WAVE)[None, :] - np.arange(360)[:, None]))


# ----------------------------------------------------------------------
# Data-file resolution 
# ----------------------------------------------------------------------
DATA_FILES = {
    "csv":     "MinorPlanetSunData.csv",
    "sp":      "corrected_sp_series.npy",
    "named":   "named_planet_angles_full.npy",
    "unnamed": "unnamed_mba_full_angles.npy",
}
_PRUNE     = {".git", "node_modules", "__pycache__", ".venv", "venv"}
_MAX_DEPTH = 4

def _candidate_dirs():
    here = Path(__file__).resolve().parent
    cwd  = Path.cwd().resolve()
    out, seen = [], set()
    for b in (here, cwd):
        for d in (b, b.parent, b.parent.parent,
                  b / "YourFolder", b.parent / "YourFolder"):
            d = d.resolve()
            if d not in seen and d.is_dir():
                seen.add(d)
                out.append(d)
    return out

def _find_data_file(filename, required=True):
    for d in _candidate_dirs():
        cand = d / filename
        if cand.is_file():
            return cand
    for start in {Path(__file__).resolve().parent, Path.cwd().resolve()}:
        for dirpath, dirnames, filenames in os.walk(start):
            rel = Path(dirpath).relative_to(start)
            if len(rel.parts) >= _MAX_DEPTH:
                dirnames[:] = []
                continue
            dirnames[:] = [d for d in dirnames
                           if d not in _PRUNE and not d.startswith(".")]
            if filename in filenames:
                return Path(dirpath) / filename
    if required:
        raise FileNotFoundError(
            f"\n  Could not find '{filename}'.\n"
            f"  Place conjunction_vtest.py and acf_significance.py in the "
            f"repo's PythonReanalysis/\n  folder, or keep the data files "
            f"within {_MAX_DEPTH - 1} folders of this script.")
    return None


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def v_statistic(wave):
    """V-test statistic for one wave evaluated at ALL 360 directions.
       Returns an array of shape (360,)."""
    total = wave.sum()
    if total <= 0:
        return np.zeros(360)
    # Dot product computes sum(w * cos) for all 360 reference directions at once
    return np.dot(_COS, wave) / total


def bf_from_p(p):
    """Pragmatic p-to-Bayes-factor conversion, 1/(2p)."""
    return 1.0 / (2.0 * p)


def bf_label(bf):
    """Verbal label for a Bayes factor on the Jeffreys scale."""
    if not np.isfinite(bf):
        return "n/a"
    if bf < 1:
        return "supports null"
    if bf < 3:
        return "anecdotal"
    if bf < 10:
        return "moderate"
    if bf < 30:
        return "strong"
    if bf < 100:
        return "very strong"
    return "extreme"


def fmt(x, w, p):
    return f"{x:>{w}.{p}f}" if np.isfinite(x) else f"{'n/a':>{w}}"


# ----------------------------------------------------------------------
def main():
    rng = np.random.default_rng(base.SEED)
    print("=" * 70)
    print("CONJUNCTION V-TEST  —  evaluating ALL 360 directions (0 to 359 deg)")
    print("=" * 70)

    # ---- locate data files ----
    print("\nLocating data files ...")
    CSV_FILE    = _find_data_file(DATA_FILES["csv"])
    SP_NPY      = _find_data_file(DATA_FILES["sp"])
    NAMED_NPY   = _find_data_file(DATA_FILES["named"])
    UNNAMED_NPY = _find_data_file(DATA_FILES["unnamed"], required=False)
    OUT_DIR     = SP_NPY.parent
    print(f"  {'CSV':<14}: {CSV_FILE}")
    print(f"  {'SP series':<14}: {SP_NPY}")
    print(f"  {'named angles':<14}: {NAMED_NPY}")
    print(f"  {'unnamed pool':<14}: "
          f"{UNNAMED_NPY if UNNAMED_NPY else '(not found - contrast skipped)'}")
    print(f"  {'output dir':<14}: {OUT_DIR}")

    # ---- data ----
    body = base.load_csv(CSV_FILE)
    articles_raw = body[:, 0::2].T
    n_planets_full, n_dates = articles_raw.shape
    named_full = np.load(NAMED_NPY)
    valid = ~np.any(np.isnan(named_full), axis=1)
    named_angles = named_full[valid]
    articles = articles_raw[valid]
    n_planets = articles.shape[0]
    sp_series = np.load(SP_NPY)
    n_sp = sp_series.shape[0]
    norm_articles = np.array(
        [base.normalize_l2(articles[i].astype(float)) for i in range(n_planets)],
        dtype=np.float64)
    print(f"\nHorizons-verified planets : {n_planets}/{n_planets_full}")
    print(f"Directions evaluated      : 360 (0 to 359 degrees)")
    print(f"Monte-Carlo replicates    : {N_MC:,} per null")

    # named-planet angle tables
    ang2d = np.empty((n_sp, n_planets, n_dates), dtype=np.int64)
    for j in range(n_sp):
        ang = (named_angles - sp_series[j][None, :]) % 360.0
        ang2d[j] = np.ceil(ang).astype(np.int64) % 360
    ang_flat = [ang2d[j].reshape(-1) for j in range(n_sp)]

    # ---- observed V ----
    obs_v = np.zeros((n_sp, 360))
    for j in range(n_sp):
        wave = base.wave_from_pairs(ang_flat[j], norm_articles.reshape(-1))
        obs_v[j, :] = v_statistic(wave)

    # ---- null distributions ----
    def run_null(kind, unnamed_ang=None):
        nv = np.zeros((N_MC, n_sp, 360))
        row = np.arange(n_planets)[:, None]
        cols = np.arange(n_dates)[None, :]
        t = time.time()
        rep = max(1, N_MC // 10)
        for i in range(N_MC):
            if i % rep == 0:
                el = time.time() - t
                eta = el / max(i, 1) * (N_MC - i)
                print(f"      iter {i:6,d}/{N_MC:,}  "
                      f"elapsed {el:4.0f}s  ETA {eta:4.0f}s")
            if kind == "compound":
                perm = rng.permutation(n_planets)
                sw = norm_articles[perm]
                sh = rng.integers(0, n_dates, n_planets)
                flat = sw[row, (cols + sh[:, None]) % n_dates].reshape(-1)
                for j in range(n_sp):
                    w = base.wave_from_pairs(ang_flat[j], flat)
                    nv[i, j, :] = v_statistic(w)
            else:  # unnamed-MBA
                perm = rng.permutation(n_planets)
                flat = norm_articles.reshape(-1)
                for j in range(n_sp):
                    w = base.wave_from_pairs(
                        unnamed_ang[j][perm].reshape(-1), flat)
                    nv[i, j, :] = v_statistic(w)
        print(f"      {kind} null done in {time.time()-t:.0f}s")
        return nv

    print(f"\n[1/2] Compound-null MC ...")
    cmp_v = run_null("compound")

    print(f"\n[2/2] Unnamed-MBA-null MC ...")
    if UNNAMED_NPY is not None:
        pool = np.load(UNNAMED_NPY)
        if pool.shape[0] >= n_planets:
            pool = pool[:n_planets]
            unnamed_ang = np.empty((n_sp, n_planets, n_dates), dtype=np.int64)
            for j in range(n_sp):
                ang = (pool - sp_series[j][None, :]) % 360.0
                unnamed_ang[j] = np.ceil(ang).astype(np.int64) % 360
            mba_v = run_null("unnamed_mba", unnamed_ang)
            have_mba = True
        else:
            have_mba = False
            print(f"      unnamed pool has only {pool.shape[0]} bodies "
                  f"(< {n_planets}); skipping unnamed-MBA contrast.")
    else:
        have_mba = False
        print(f"      {DATA_FILES['unnamed']} not found "
              f"-- skipping the unnamed-MBA contrast.")

    # ---- statistics ----
    def stats(obs, null):
        """Per-SP and omnibus statistics for the V-statistic under one null.
           obs: (n_sp, 360), null: (N_MC, n_sp, 360)."""
        mu = null.mean(axis=0)                                 # (n_sp, 360)
        sd = null.std(axis=0, ddof=1)                          # (n_sp, 360)
        d  = (obs - mu) / sd                                   # signed effect
        dev_obs  = np.abs(obs - mu)
        dev_null = np.abs(null - mu)
        p  = (1 + (dev_null >= dev_obs).sum(axis=0)) / (1 + N_MC)   # two-sided
        
        # Omnibus across the 11 special points, per direction
        T_obs  = ((obs - mu) ** 2 / sd ** 2).sum(axis=0)       # (360,)
        T_null = ((null - mu) ** 2 / sd ** 2).sum(axis=1)      # (N_MC, 360)
        p_omni = (1 + (T_null >= T_obs).sum(axis=0)) / (1 + N_MC)
        
        return dict(mu=mu, sd=sd, d=d, p=p,
                    T_obs=T_obs, T_null_mean=T_null.mean(axis=0),
                    p_omni=p_omni)

    sc = stats(obs_v, cmp_v)
    sm = stats(obs_v, mba_v) if have_mba else None

    # ---- console preview ----
    deg = 0  # We will just preview 0 degrees in the terminal
    print("\n" + "=" * 78)
    print(f"PREVIEW: V-STATISTIC FOR 0 DEG (CONJUNCTION)")
    print("  Full results for all 0-359 deg are saved in the CSVs.")
    print("  d = (V - null mean) / null SD ;  p two-sided ;  BF = 1/(2p)")
    print("=" * 78)
    print(f"  {'SP':<12}{'V_obs':>11}"
          f"{'d(cmp)':>9}{'p(cmp)':>9}{'BF(cmp)':>9}"
          f"{'d(mba)':>9}{'p(mba)':>9}{'BF(mba)':>9}")
    print("  " + "-" * 74)
    
    for j in range(n_sp):
        dc, pc = sc["d"][j, deg], sc["p"][j, deg]
        bc = bf_from_p(pc)
        if sm is not None:
            dm, pm = sm["d"][j, deg], sm["p"][j, deg]
            bm = bf_from_p(pm)
        else:
            dm = pm = bm = np.nan
            
        print(f"  {LABELS[j]:<12}{obs_v[j, deg]:>11.5f}"
              f"{fmt(dc,9,2)}{fmt(pc,9,4)}{fmt(bc,9,1)}"
              f"{fmt(dm,9,2)}{fmt(pm,9,4)}{fmt(bm,9,1)}")
              
    print("  " + "-" * 74)
    bc_o = bf_from_p(sc["p_omni"][deg])
    print(f"  Omnibus (\u03a3 d\u00b2)  compound: T = {sc['T_obs'][deg]:.1f}  "
          f"(null mean {sc['T_null_mean'][deg]:.1f})   "
          f"p = {sc['p_omni'][deg]:.4f}   BF = {bc_o:.1f}  [{bf_label(bc_o)}]")
    if sm is not None:
        bm_o = bf_from_p(sm["p_omni"][deg])
        print(f"  Omnibus (\u03a3 d\u00b2)  unnamed : T = {sm['T_obs'][deg]:.1f}  "
              f"(null mean {sm['T_null_mean'][deg]:.1f})   "
              f"p = {sm['p_omni'][deg]:.4f}   BF = {bm_o:.1f}  [{bf_label(bm_o)}]")

    # ---- generate full CSV data for all 360 degrees ----
    csv_per_sp = []
    csv_omni = []
    
    for d in range(360):
        # Per SP rows
        for j in range(n_sp):
            bc = bf_from_p(sc["p"][j, d])
            bm = bf_from_p(sm["p"][j, d]) if sm is not None else np.nan
            csv_per_sp.append([
                LABELS[j], d, f"{obs_v[j, d]:.8f}",
                f"{sc['mu'][j, d]:.8f}", f"{sc['sd'][j, d]:.8f}",
                f"{sc['d'][j, d]:.6f}", f"{sc['p'][j, d]:.6f}", f"{bc:.6f}",
                (f"{sm['mu'][j, d]:.8f}" if sm else ""),
                (f"{sm['sd'][j, d]:.8f}" if sm else ""),
                (f"{sm['d'][j, d]:.6f}" if sm else ""),
                (f"{sm['p'][j, d]:.6f}" if sm else ""),
                (f"{bm:.6f}" if sm else "")])
                
        # Omnibus rows
        bc_o = bf_from_p(sc["p_omni"][d])
        bm_o = bf_from_p(sm["p_omni"][d]) if sm is not None else np.nan
        csv_omni.append([
            d,
            f"{sc['T_obs'][d]:.6f}", f"{sc['T_null_mean'][d]:.6f}",
            f"{sc['p_omni'][d]:.6f}", f"{bc_o:.6f}",
            (f"{sm['T_obs'][d]:.6f}" if sm else ""),
            (f"{sm['T_null_mean'][d]:.6f}" if sm else ""),
            (f"{sm['p_omni'][d]:.6f}" if sm else ""),
            (f"{bm_o:.6f}" if sm else "")])

    # ---- write CSVs ----
    out1 = OUT_DIR / "conjunction_vtest_per_sp.csv"
    with open(out1, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["special_point", "direction_deg", "v_obs",
                    "cmp_null_mean", "cmp_null_sd", "cmp_d", "cmp_p", "cmp_bf",
                    "mba_null_mean", "mba_null_sd", "mba_d", "mba_p", "mba_bf"])
        w.writerows(csv_per_sp)

    out2 = OUT_DIR / "conjunction_vtest_omnibus.csv"
    with open(out2, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["direction_deg",
                    "cmp_T_obs", "cmp_T_null_mean", "cmp_p", "cmp_bf",
                    "mba_T_obs", "mba_T_null_mean", "mba_p", "mba_bf"])
        w.writerows(csv_omni)

    print(f"\nSaved all 360 directions -> {out1}")
    print(f"Saved all 360 directions -> {out2}")

if __name__ == "__main__":
    main()
