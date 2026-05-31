"""
all_sps_vs_gauquelin_extracted.py
   — Same protocol as mars_vs_gauquelin_extracted_shape.py, but applied
     to every non-Mesarthim special point.

  Workflow
  --------
  1. Extract the Gauquelin reference shape from
        ~/Downloads/gauquelinzonesrotated.gif
     by tracing radial distances from the cross-centre at 360 integer
     angles in (0°=top, CW) convention.  Mean-centre and peak-normalise.
  2. For each SP j ∈ {Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn,
     Uranus, Neptune, Pluto}:
        a. Build wave_j[k] (bin k = angle to SP_j).
        b. Rotate by φ̄_j (circular mean of SP_j_lon) into Mes-frame.
        c. Gaussian-smooth (σ = 30 bin, mode='wrap'), mean-centre.
        d. Take residual_j = SP_j_smooth_centred − Mes_smooth_centred.
        e. Peak-normalise to ±1.
        f. Pearson r vs extracted Gauquelin shape  (and best circular
           shift to gauge alignment robustness).
  3. Render a 2×5 polar grid of overlays + a ranking table + summary
     plots.
"""

import csv as csv_mod
import warnings
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

warnings.filterwarnings("ignore")

DATA_DIR  = Path(__file__).parent
OUT_DIR   = Path(__file__).parent
CSV_FILE  = DATA_DIR / "MinorPlanetSunData.csv"
SP_NPY    = OUT_DIR / "corrected_sp_series.npy"
NAMED_NPY = OUT_DIR / "named_planet_angles_full.npy"
GAUQ_GIF  = Path(__file__).parent / "gauquelinzonesrotated.gif"

N_WAVE     = 360
SIGMA      = 30
LABELS = [
    "Mesarthim", "Sun", "Moon", "Mercury", "Venus",
    "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Rahu",
]
MESARTHIM_IDX = 0
SP_COLORS = [
    "#666666", "#e6a700", "#888aa8", "#7a5cb1", "#d96aa6",
    "#c0392b", "#b07a1a", "#3a7a52", "#2c8aa1", "#26426b", "#5a2c80", "#6a1b9a",
]
ASTRO_LAGS = [24, 30, 36, 40, 45, 51, 60, 72, 90, 120, 180]


def _parse(s):
    try:    return float(s.strip())
    except ValueError: return float("nan")


def load_csv(path):
    with open(path, newline="") as f:
        rows = list(csv_mod.reader(f))[:-1]
    return np.array([[_parse(c) for c in r[1:]] for r in rows[1:]], dtype=float)


def normalize_l2(v):
    n = np.linalg.norm(v); return v / n if n > 0 else v.copy()


def wave_from_pairs(angle_flat, value_flat, n_deg=N_WAVE):
    sums = np.bincount(angle_flat, weights=value_flat, minlength=n_deg)
    cnts = np.bincount(angle_flat,                     minlength=n_deg).clip(1)
    return sums / cnts


def circular_mean_deg(a):
    return np.degrees(np.angle(np.mean(np.exp(1j * np.radians(a))))) % 360.0


def extract_gauquelin_shape(gif_path):
    im = Image.open(str(gif_path)).convert("L")
    arr = np.array(im).astype(np.int16)
    H, W = arr.shape
    black = arr < 60
    gray  = (arr >= 60) & (arr < 200)
    cx = int(np.argmax(gray.sum(axis=0)))
    cy = int(np.argmax(gray.sum(axis=1)))
    max_r = int(0.55 * min(H, W))
    radii = np.zeros(N_WAVE, dtype=np.float64)
    for theta_deg in range(N_WAVE):
        rad = np.radians(theta_deg)
        dx = np.sin(rad); dy = -np.cos(rad)
        last_black_r = 0.0
        for rr in range(1, max_r):
            x = int(round(cx + dx * rr))
            y = int(round(cy + dy * rr))
            if x < 0 or x >= W or y < 0 or y >= H: break
            if black[y, x]: last_black_r = float(rr)
        radii[theta_deg] = last_black_r
    zero_mask = radii == 0
    if zero_mask.any():
        good_idx = np.where(~zero_mask)[0]
        bad_idx  = np.where(zero_mask)[0]
        radii[zero_mask] = np.interp(bad_idx, good_idx, radii[good_idx], period=N_WAVE)
    # circular median filter
    radii = median_filter(np.concatenate([radii[-7:], radii, radii[:7]]), size=5)[7:-7]
    return radii, (cx, cy)


def main():
    print(f"Loading Gauquelin reference: {GAUQ_GIF}")
    gauq_radii, centre = extract_gauquelin_shape(GAUQ_GIF)
    gauq_norm = (gauq_radii - gauq_radii.mean())
    gauq_norm = gauq_norm / np.max(np.abs(gauq_norm))
    pk_g = int(np.argmax(gauq_norm)); tr_g = int(np.argmin(gauq_norm))
    print(f"  Extracted: peak @ {pk_g}°  ·  trough @ {tr_g}°  ·  centre {centre}")

    body = load_csv(CSV_FILE)
    articles_raw = body[:, 0::2].T
    named_angles_full = np.load(NAMED_NPY)
    valid_mask = ~np.any(np.isnan(named_angles_full), axis=1)
    n_planets = int(valid_mask.sum())
    named_angles = named_angles_full[valid_mask]
    articles = articles_raw[valid_mask]
    sp_series = np.load(SP_NPY)
    n_sp = sp_series.shape[0]
    sp_lon_t = (sp_series - sp_series[MESARTHIM_IDX:MESARTHIM_IDX + 1]) % 360.0
    phi_bar = np.array([circular_mean_deg(sp_lon_t[j]) for j in range(n_sp)])
    norm_articles = np.array(
        [normalize_l2(articles[i].astype(float)) for i in range(n_planets)],
        dtype=np.float32
    )

    print(f"\nBuilding per-SP waves (n={n_planets} Horizons-verified planets)…")
    waves_sp = np.zeros((n_sp, N_WAVE))
    for j in range(n_sp):
        sp = sp_series[j]
        ang = (named_angles - sp[None, :]) % 360.0
        ang_int = (np.ceil(ang).astype(np.int32) % 360).reshape(-1).astype(np.int64)
        waves_sp[j] = wave_from_pairs(ang_int, norm_articles.reshape(-1))
    waves_mes_frame = np.zeros_like(waves_sp)
    for j in range(n_sp):
        waves_mes_frame[j] = np.roll(waves_sp[j], int(np.round(phi_bar[j])) % N_WAVE)
    waves_smooth = np.array([gaussian_filter1d(w, SIGMA, mode="wrap")
                              for w in waves_mes_frame])
    waves_centered = waves_smooth - waves_smooth.mean(axis=1, keepdims=True)
    mes_curve = waves_centered[MESARTHIM_IDX]
    residuals = waves_centered - mes_curve[None, :]   # Mes residual = 0 by construction

    # Peak-normalise each non-Mesarthim residual
    sp_indices = [j for j in range(n_sp) if j != MESARTHIM_IDX]
    norm_residuals = {}
    summary = []
    for j in sp_indices:
        r = residuals[j]
        peak_amp = float(np.max(np.abs(r)))
        rn = r / peak_amp if peak_amp > 0 else r
        norm_residuals[j] = rn
        pk_m = int(np.argmax(rn)); tr_m = int(np.argmin(rn))
        pearson_r = float(np.corrcoef(rn, gauq_norm)[0, 1])
        rs = np.array([np.corrcoef(np.roll(rn, sh), gauq_norm)[0, 1]
                        for sh in range(N_WAVE)])
        best_shift = int(np.argmax(rs)); best_r = float(rs[best_shift])
        summary.append({
            "j": j,
            "label": LABELS[j],
            "peak_bin": pk_m,
            "trough_bin": tr_m,
            "peak_amp_raw": peak_amp,
            "r_pearson": pearson_r,
            "r_pearson_sq": pearson_r ** 2,
            "best_shift": best_shift,
            "r_best_shift": best_r,
        })

    # Ranking
    summary.sort(key=lambda x: -x["r_pearson"])
    print(f"\n══════════ RANKING: Pearson r(SP residual σ={SIGMA}, extracted Gauquelin) ══════════")
    print(f"  {'rank':<5} {'SP':<10} {'peak θ':>7} {'trough θ':>9} "
          f"{'r':>9} {'R²':>8} {'best_shift':>11} {'r_best':>8}")
    print("  " + "-" * 75)
    for rank, s in enumerate(summary, start=1):
        # convert best_shift to ±180° signed
        sh = s["best_shift"]; sh_s = sh - 360 if sh > 180 else sh
        print(f"  {rank:<5} {s['label']:<10} {s['peak_bin']:>6d}° {s['trough_bin']:>8d}° "
              f"{s['r_pearson']:>+9.4f} {s['r_pearson_sq']:>8.4f} "
              f"{sh_s:>+10d}° {s['r_best_shift']:>+8.4f}")

    # ── PLOTTING: 2×5 polar grid of overlays ─────────────────── #
    theta = np.radians(np.arange(N_WAVE))
    theta_closed = np.concatenate([theta, theta[:1]])
    rg = 0.85 + 0.45 * gauq_norm
    rg_closed = np.concatenate([rg, rg[:1]])

    fig = plt.figure(figsize=(28, 13))
    # plot in ranking order (so the strongest matches come first)
    for idx, s in enumerate(summary):
        ax = fig.add_subplot(3, 4, idx + 1, projection="polar")
        rn = norm_residuals[s["j"]]
        rd = 0.85 + 0.45 * rn
        rd_closed = np.concatenate([rd, rd[:1]])
        # Background: reference circle, spokes, crosshairs
        ax.plot(theta_closed, np.ones_like(theta_closed) * 0.85,
                color="gray", lw=1.0, alpha=0.7)
        for ang_deg in range(0, 360, 45):
            ax.plot([np.radians(ang_deg), np.radians(ang_deg)], [0, 1.4],
                    color="gray", lw=0.5, alpha=0.5)
        for ang_deg in [0, 90, 180, 270]:
            ax.plot([np.radians(ang_deg), np.radians(ang_deg)], [0, 1.4],
                    color="dimgray", lw=1.2, alpha=0.7)
        # Gauquelin shape (faint orange)
        ax.plot(theta_closed, rg_closed, color="#d97a00", lw=2.0,
                alpha=0.55, label="Gauquelin (extracted)")
        ax.fill(theta_closed, rg_closed, color="#fbd29a", alpha=0.20)
        # SP residual (bold, SP-coloured)
        ax.plot(theta_closed, rd_closed, color=SP_COLORS[s["j"]],
                lw=2.6, zorder=10, label=f"{s['label']} (σ={SIGMA})")
        # Peak/trough markers
        ax.scatter([np.radians(s["peak_bin"])], [rd[s["peak_bin"]]],
                    color="darkred", s=50, zorder=11,
                    edgecolors="white", linewidths=0.8)
        ax.scatter([np.radians(s["trough_bin"])], [rd[s["trough_bin"]]],
                    color="darkblue", s=50, zorder=11,
                    edgecolors="white", linewidths=0.8)
        # Peak of Gauquelin for reference
        ax.scatter([np.radians(pk_g)], [rg[pk_g]],
                    marker="*", color="#d97a00", s=110, zorder=12,
                    edgecolors="white", linewidths=0.6)

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(np.radians(np.arange(0, 360, 45)))
        ax.set_xticklabels(["0", "45", "90", "135", "180", "225", "270", "315"],
                            fontsize=7)
        ax.set_ylim(0, 1.4)
        ax.set_yticklabels([])
        ax.set_title(
            f"#{idx+1}  {s['label']}\n"
            f"peak {s['peak_bin']}° ·  r = {s['r_pearson']:+.3f}\n"
            f"best shift {s['best_shift'] - 360 if s['best_shift'] > 180 else s['best_shift']:+d}°  →  r = {s['r_best_shift']:+.3f}",
            fontsize=10, pad=10
        )

    fig.suptitle(
        f"All non-Mesarthim SPs vs EXTRACTED Gauquelin shape\n"
        f"Each panel: SP − Mes residual (σ={SIGMA}, normalised) overlaid on extracted Gauquelin (faint orange).\n"
        f"Panels ordered by Pearson r (best match first).  Orange star = Gauquelin peak @ {pk_g}°.  "
        f"n = {n_planets} Horizons-verified named planets.",
        fontsize=12, y=1.005
    )
    plt.tight_layout()
    out = OUT_DIR / "all_sps_vs_gauquelin_extracted.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot → {out}")

    # ── Summary bar chart of Pearson r per SP ──────────────────── #
    fig2, ax2 = plt.subplots(figsize=(13, 6))
    names    = [s["label"] for s in summary]
    rvals    = [s["r_pearson"] for s in summary]
    rbestvals = [s["r_best_shift"] for s in summary]
    cols     = [SP_COLORS[s["j"]] for s in summary]
    xs = np.arange(len(names))
    bars = ax2.bar(xs - 0.18, rvals, width=0.36, color=cols, edgecolor="black",
                    linewidth=0.8, label="Pearson r (zero shift)")
    bars2 = ax2.bar(xs + 0.18, rbestvals, width=0.36, color=cols, edgecolor="black",
                     linewidth=0.8, alpha=0.5, label="Pearson r (best shift)")
    for x, r, rb in zip(xs, rvals, rbestvals):
        ax2.text(x - 0.18, r + (0.02 if r >= 0 else -0.04),
                  f"{r:+.3f}", ha="center", fontsize=8)
        ax2.text(x + 0.18, rb + (0.02 if rb >= 0 else -0.04),
                  f"{rb:+.3f}", ha="center", fontsize=8, alpha=0.7)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(names, rotation=0, fontsize=10)
    ax2.set_ylim(-0.8, 1.0)
    ax2.axhline(0, color="black", lw=0.6)
    ax2.set_ylabel("Pearson r vs extracted Gauquelin shape", fontsize=11)
    ax2.set_title(f"Match between each SP's σ={SIGMA} residual and the extracted Gauquelin curve\n"
                   f"(ranking left → right by Pearson r at zero shift)",
                   fontsize=11)
    ax2.grid(alpha=0.3, axis="y")
    ax2.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    out2 = OUT_DIR / "all_sps_vs_gauquelin_rank.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot → {out2}")

    # ── Linear overlay panel: each SP residual + Gauquelin (stacked) ── #
    fig3, axes = plt.subplots(12, 1, figsize=(16, 26), sharex=True)
    deg = np.arange(N_WAVE)
    for ax, s in zip(axes, summary):
        rn = norm_residuals[s["j"]]
        ax.plot(deg, gauq_norm, color="#d97a00", lw=1.4, alpha=0.55,
                 label="Gauquelin (extracted)")
        ax.fill_between(deg, 0, gauq_norm, where=(gauq_norm > 0),
                         color="#fbd29a", alpha=0.30)
        ax.plot(deg, rn, color=SP_COLORS[s["j"]], lw=1.8,
                 label=f"{s['label']} (σ={SIGMA})")
        ax.axhline(0, color="gray", lw=0.4)
        for lag in ASTRO_LAGS:
            ax.axvline(lag, color="gray", lw=0.2, alpha=0.4)
        ax.set_xlim(0, 360)
        ax.set_ylim(-1.2, 1.2)
        ax.set_ylabel(f"{s['label']}\n"
                       f"r = {s['r_pearson']:+.3f}", fontsize=9)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Bin α (degrees from Mesarthim)", fontsize=10)
    fig3.suptitle("Linear overlay — each SP residual (σ={SIGMA}, normalised) vs extracted Gauquelin"
                   .format(SIGMA=SIGMA),
                   fontsize=11, y=1.001)
    plt.tight_layout()
    out3 = OUT_DIR / "all_sps_vs_gauquelin_linear.png"
    plt.savefig(out3, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot → {out3}")

    # ── CSV ─────────────────────────────────────────────────────── #
    with open(OUT_DIR / "all_sps_vs_gauquelin_extracted_summary.csv", "w", newline="") as f:
        w = csv_mod.writer(f)
        w.writerow(["rank", "SP", "peak_bin_deg", "trough_bin_deg",
                    "raw_peak_amp", "pearson_r_zero_shift", "R_sq",
                    "best_circular_shift_deg", "pearson_r_at_best_shift"])
        for rank, s in enumerate(summary, 1):
            sh = s["best_shift"]; sh_s = sh - 360 if sh > 180 else sh
            w.writerow([rank, s["label"], s["peak_bin"], s["trough_bin"],
                        f"{s['peak_amp_raw']:.6e}",
                        f"{s['r_pearson']:+.6f}", f"{s['r_pearson_sq']:.6f}",
                        sh_s, f"{s['r_best_shift']:+.6f}"])
    print(f"  CSV  → {OUT_DIR / 'all_sps_vs_gauquelin_extracted_summary.csv'}")

    # Full per-bin residuals
    with open(OUT_DIR / "all_sps_vs_gauquelin_extracted_curves.csv", "w", newline="") as f:
        w = csv_mod.writer(f)
        header = ["bin_deg", "gauquelin_extracted_norm"] + \
                 [f"{LABELS[j]}_residual_sigma{SIGMA}_norm" for j in sp_indices]
        w.writerow(header)
        for k in range(N_WAVE):
            row = [k, f"{gauq_norm[k]:+.6e}"] + \
                  [f"{norm_residuals[j][k]:+.6e}" for j in sp_indices]
            w.writerow(row)
    print(f"  CSV  → {OUT_DIR / 'all_sps_vs_gauquelin_extracted_curves.csv'}")


if __name__ == "__main__":
    main()
