"""
acf_plots_all_sps.py
   — Per-SP autocorrelation-function plots over all lags 0..359°, with
     the white-noise 95% confidence band shaded and the astrological
     periods (24°, 30°, 36°, 40°, 45°, 51°, 60°, 72°, 90°, 120°, 180°)
     marked pictorially.

  Methodology
  -----------
  • Each per-SP wave wave_j[k]  is built with bin k = angle of planet
    to SP_j, summed over all (planet, day) pairs, normalised by count.
  • Biased ACF is computed on the wave:
        ACF[ℓ] = Σ_i (x[i]·x[(i+ℓ) mod N]) / Σ_i x[i]²    where x = wave − mean
    (matches Mathematica CorrelationFunction and the project's other
    ACF code).
  • The 95% confidence band for circular ACF under white-noise H₀:
        ± 1.96 / √N   with N = 360.
  • Astrological aspects are drawn as vertical spokes with labels and a
    coloured dot whose colour reflects whether |ACF(ℓ)| > the CI bound.
  • ACF is invariant under circular rotation, so the SP-frame and the
    Mes-frame waves give identical ACFs — we use the SP-frame wave (the
    canonical version) here.

  Outputs
  -------
  1. 4×3 grid of linear ACF plots, one panel per SP  (acf_grid_linear.png).
  2. 4×3 grid of polar ACF plots, one panel per SP   (acf_grid_polar.png).
  3. Per-SP standalone ACF plots (acf_<sp>.png) for high-resolution
     individual viewing.
  4. CSV of ACF at all 360 lags for all 11 SPs.
"""

import csv as csv_mod
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

DATA_DIR  = Path(__file__).parent
OUT_DIR   = Path(__file__).parent
CSV_FILE  = DATA_DIR / "MinorPlanetSunData.csv"
SP_NPY    = OUT_DIR / "corrected_sp_series.npy"
NAMED_NPY = OUT_DIR / "named_planet_angles_full.npy"

N_WAVE = 360
ACF_THRESHOLD = 1.96 / np.sqrt(N_WAVE)

LABELS = [
    "Mesarthim", "Sun", "Moon", "Mercury", "Venus",
    "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Rahu",
]
MESARTHIM_IDX = 0

ASTRO_LAGS = [24, 30, 36, 40, 45, 51, 60, 72, 90, 120, 180]
ASTRO_NAMES = {
    24:  "quindec.",
    30:  "semi-sxt",
    36:  "decile",
    40:  "novile",
    45:  "semi-sq",
    51:  "septile",
    60:  "sextile",
    72:  "quintile",
    90:  "square",
    120: "trine",
    180: "opp.",
}

SP_COLORS = [
    "#666666", "#e6a700", "#888aa8", "#7a5cb1", "#d96aa6",
    "#c0392b", "#b07a1a", "#3a7a52", "#2c8aa1", "#26426b", "#5a2c80", "#6a1b9a",
]


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


def biased_acf_circular(wave):
    """Biased circular ACF at all lags 0..N-1.  Returns shape (N,)."""
    x = np.asarray(wave, dtype=np.float64) - wave.mean()
    n = len(x)
    var = (x * x).sum()
    out = np.empty(n)
    if var == 0.0:
        return out
    # Use FFT for efficiency: ACF = IFFT(|FFT(x)|²) / var  (circular)
    F = np.fft.fft(x)
    out = np.real(np.fft.ifft(F * np.conj(F))) / var
    return out


def build_waves(named_angles, sp_series, norm_articles, n_sp):
    waves = np.zeros((n_sp, N_WAVE))
    for j in range(n_sp):
        sp = sp_series[j]
        ang = (named_angles - sp[None, :]) % 360.0
        ang_int = (np.ceil(ang).astype(np.int32) % 360).reshape(-1).astype(np.int64)
        waves[j] = wave_from_pairs(ang_int, norm_articles.reshape(-1))
    return waves


def draw_linear_panel(ax, acf, sp_name, color, show_xticks=True, show_legend=False):
    deg = np.arange(N_WAVE)
    # Significance shading: anywhere |ACF| > CI is highlighted
    ax.fill_between(deg, -ACF_THRESHOLD, ACF_THRESHOLD,
                    color="lightgray", alpha=0.35, label="95% CI band")
    ax.axhline(0, color="black", lw=0.5)
    ax.axhline(ACF_THRESHOLD, color="gray", linestyle="--", lw=0.6, alpha=0.7)
    ax.axhline(-ACF_THRESHOLD, color="gray", linestyle="--", lw=0.6, alpha=0.7)
    # ACF curve
    ax.plot(deg, acf, color=color, lw=1.4, alpha=0.95, label=f"ACF ({sp_name})")
    # Shade ACF > CI in red (positive significant) and ACF < -CI in blue (negative significant)
    above = acf > ACF_THRESHOLD
    below = acf < -ACF_THRESHOLD
    ax.fill_between(deg, ACF_THRESHOLD, acf, where=above,
                    color=color, alpha=0.30, interpolate=True)
    ax.fill_between(deg, -ACF_THRESHOLD, acf, where=below,
                    color="firebrick", alpha=0.30, interpolate=True)
    # Astrological lag spokes
    y_top = max(0.55, max(acf[1:].max() * 1.1, ACF_THRESHOLD * 2))
    y_bot = min(-0.55, min(acf[1:].min() * 1.1, -ACF_THRESHOLD * 2))
    for lag in ASTRO_LAGS:
        v = acf[lag]
        is_sig = abs(v) > ACF_THRESHOLD
        spoke_color = "darkred" if v > 0 else "darkblue"
        spoke_alpha = 0.9 if is_sig else 0.4
        # Spoke
        ax.plot([lag, lag], [0, v], color=spoke_color,
                lw=1.6 if is_sig else 0.7,
                alpha=spoke_alpha, zorder=4)
        # Dot
        ax.scatter([lag], [v], color=spoke_color,
                   s=70 if is_sig else 28, zorder=6,
                   edgecolors="white", linewidths=0.8 if is_sig else 0.4)
        # Label above/below the dot
        label_y = v + 0.04 * np.sign(v) if v != 0 else 0.04
        ax.text(lag, label_y, f"{lag}°\n{ASTRO_NAMES[lag]}\n{v:+.3f}",
                 ha="center", va="bottom" if v > 0 else "top",
                 fontsize=6.5, color=spoke_color,
                 fontweight="bold" if is_sig else "normal",
                 alpha=0.95 if is_sig else 0.7)
    ax.set_xlim(0, 359)
    # Y-limits: tight to data, but at least ±CI*2
    ymax_data = max(abs(acf[1:].max()), abs(acf[1:].min()))
    yextent = max(ymax_data * 1.30, ACF_THRESHOLD * 3.0)
    ax.set_ylim(-yextent, yextent)
    ax.set_ylabel(f"{sp_name}\nACF", fontsize=10, color=color, fontweight="bold")
    ax.grid(alpha=0.3)
    if show_xticks:
        ax.set_xticks(list(range(0, 361, 30)))
        ax.set_xticklabels([f"{a}°" for a in range(0, 361, 30)], fontsize=8)
    else:
        ax.set_xticks(list(range(0, 361, 30)))
        ax.set_xticklabels([])
    if show_legend:
        # Compact legend
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        handles = [
            Patch(facecolor="lightgray", alpha=0.5,
                   label=f"95% CI band  (±{ACF_THRESHOLD:.4f})"),
            Line2D([0],[0], color="gray", lw=1, linestyle="--",
                   label="CI threshold"),
            Line2D([0],[0], marker="o", color="darkred", lw=0,
                   markersize=8, label="Astro lag (ACF > +CI)"),
            Line2D([0],[0], marker="o", color="darkblue", lw=0,
                   markersize=8, label="Astro lag (ACF < −CI)"),
            Line2D([0],[0], marker="o", color="dimgray", lw=0,
                   markersize=5, label="Astro lag (within CI)"),
        ]
        ax.legend(handles=handles, loc="upper right", fontsize=7,
                  framealpha=0.85)


def draw_polar_panel(ax, acf, sp_name, color):
    theta = np.radians(np.arange(N_WAVE))
    theta_closed = np.concatenate([theta, theta[:1]])
    r = acf.copy()
    r_closed = np.concatenate([r, r[:1]])
    # 95% CI band as a faint ring
    ax.fill_between(theta_closed, -ACF_THRESHOLD, ACF_THRESHOLD,
                    color="lightgray", alpha=0.35)
    ax.plot(theta_closed, np.ones_like(theta_closed) * ACF_THRESHOLD,
            color="gray", lw=0.6, linestyle="--", alpha=0.7)
    ax.plot(theta_closed, np.ones_like(theta_closed) * -ACF_THRESHOLD,
            color="gray", lw=0.6, linestyle="--", alpha=0.7)
    ax.plot(theta_closed, np.zeros_like(theta_closed),
            color="black", lw=0.5, alpha=0.5)
    # ACF curve
    ax.plot(theta_closed, r_closed, color=color, lw=1.5, alpha=0.95)
    # Fill positive in colour, negative in red
    ax.fill_between(theta, 0, np.where(r > 0, r, 0), color=color, alpha=0.35)
    ax.fill_between(theta, 0, np.where(r < 0, r, 0), color="firebrick", alpha=0.30)
    # Astro spokes
    for lag in ASTRO_LAGS:
        ang = np.radians(lag); v = acf[lag]
        is_sig = abs(v) > ACF_THRESHOLD
        sc = "darkred" if v > 0 else "darkblue"
        ax.plot([ang, ang], [0, v], color=sc,
                lw=1.5 if is_sig else 0.7,
                alpha=0.9 if is_sig else 0.5)
        ax.scatter([ang], [v], color=sc,
                   s=55 if is_sig else 22, zorder=6,
                   edgecolors="white", linewidths=0.6)
        # outer label
        label_r = max(abs(r.max()), 0.20) * 1.20
        ax.text(ang, label_r, f"{lag}°", ha="center", va="center",
                 fontsize=6.5,
                 color=sc, fontweight="bold" if is_sig else "normal")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_xticks(np.radians(np.arange(0, 360, 45)))
    ax.set_xticklabels(["0°","45°","90°","135°","180°","225°","270°","315°"],
                       fontsize=7)
    rmax = max(abs(r.max()), abs(r.min())) * 1.35
    ax.set_ylim(-rmax, rmax)
    ax.set_yticklabels([])
    sig_count = sum(1 for lag in ASTRO_LAGS if abs(acf[lag]) > ACF_THRESHOLD)
    ax.set_title(f"{sp_name}\nmax|ACF| = {max(abs(acf[1:])):.3f}  ·  "
                  f"{sig_count}/{len(ASTRO_LAGS)} astro lags |ACF|>CI",
                  fontsize=10, color=color, pad=12)
    ax.grid(alpha=0.25)


def main():
    body = load_csv(CSV_FILE)
    articles_raw = body[:, 0::2].T
    named_angles_full = np.load(NAMED_NPY)
    valid_mask = ~np.any(np.isnan(named_angles_full), axis=1)
    n_planets = int(valid_mask.sum())
    named_angles = named_angles_full[valid_mask]
    articles = articles_raw[valid_mask]
    sp_series = np.load(SP_NPY)
    n_sp = sp_series.shape[0]

    norm_articles = np.array(
        [normalize_l2(articles[i].astype(float)) for i in range(n_planets)],
        dtype=np.float32
    )

    print(f"Building per-SP waves (n={n_planets} Horizons-verified planets)…")
    waves = build_waves(named_angles, sp_series, norm_articles, n_sp)

    # Compute ACFs
    print(f"\nComputing biased circular ACFs for {n_sp} SPs (all lags 0..359)…")
    acfs = np.array([biased_acf_circular(waves[j]) for j in range(n_sp)])

    # Per-SP summary
    print(f"\nPer-SP ACF summary:")
    print(f"  {'SP':<12} {'max|ACF|':>10} {'@lag':>5}  "
          f"{'astro_lag_sig_count':>20}")
    for j in range(n_sp):
        ml = int(np.argmax(np.abs(acfs[j, 1:]))) + 1
        sig = sum(1 for lag in ASTRO_LAGS if abs(acfs[j, lag]) > ACF_THRESHOLD)
        print(f"  {LABELS[j]:<12} {max(abs(acfs[j, 1:])):>9.4f}  {ml:>4d}°  "
              f"{sig}/{len(ASTRO_LAGS)}")

    print(f"\n95% CI threshold (white-noise H₀): ±{ACF_THRESHOLD:.4f}")

    # ── 4×3 GRID OF LINEAR ACF PLOTS ─────────────────────────────── #
    fig, axes = plt.subplots(4, 3, figsize=(20, 22), sharex=False)
    axes_flat = axes.flatten()
    for j in range(n_sp):
        ax = axes_flat[j]
        is_last_row = j >= 9    # Mesarthim..Pluto: rows 0,1,2,3; last row = indices 9,10
        draw_linear_panel(ax, acfs[j], LABELS[j], SP_COLORS[j],
                           show_xticks=(j in (9, 10) or j == 8),
                           show_legend=(j == 0))
        if j in (9, 10) or j == 8:
            ax.set_xlabel("Lag (degrees of angular separation)", fontsize=10)
    # Hide unused cell (the 12th)
    if n_sp < 12:
        axes_flat[11].axis("off")

    fig.suptitle(
        f"Autocorrelation Function (ACF) per Special Point — all lags 0°..359°\n"
        f"Biased circular ACF on per-SP wave (bin = angle of planet from SP). "
        f"Gray shading = 95% CI band (±{ACF_THRESHOLD:.4f}, white-noise H₀).\n"
        f"Astrological aspects marked as spokes (red = positive sig., blue = negative sig.).  "
        f"n = {n_planets} Horizons-verified named planets.",
        fontsize=12, y=1.005
    )
    plt.tight_layout()
    out = OUT_DIR / "acf_grid_linear.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot → {out}")

    # ── 4×3 GRID OF POLAR ACF PLOTS ──────────────────────────────── #
    fig2 = plt.figure(figsize=(20, 22))
    for j in range(n_sp):
        ax = fig2.add_subplot(4, 3, j + 1, projection="polar")
        draw_polar_panel(ax, acfs[j], LABELS[j], SP_COLORS[j])

    fig2.suptitle(
        f"Polar ACF plots per Special Point — radial = ACF(lag), angular = lag in degrees\n"
        f"Gray ring = 95% CI band (±{ACF_THRESHOLD:.4f}).  Astrological lags marked.  "
        f"0° at top, clockwise.  n = {n_planets} planets.",
        fontsize=12, y=1.01
    )
    plt.tight_layout()
    out2 = OUT_DIR / "acf_grid_polar.png"
    plt.savefig(out2, dpi=140, bbox_inches="tight")
    plt.close()
    print(f"  Plot → {out2}")

    # ── PER-SP STANDALONE PLOTS (linear + polar combined) ─────── #
    for j in range(n_sp):
        fig3 = plt.figure(figsize=(18, 7))
        ax_lin = fig3.add_subplot(1, 2, 1)
        ax_pol = fig3.add_subplot(1, 2, 2, projection="polar")
        draw_linear_panel(ax_lin, acfs[j], LABELS[j], SP_COLORS[j],
                           show_xticks=True, show_legend=True)
        ax_lin.set_xlabel("Lag (degrees of angular separation)", fontsize=11)
        ax_lin.set_title(f"{LABELS[j]} — linear ACF", fontsize=12)
        draw_polar_panel(ax_pol, acfs[j], LABELS[j], SP_COLORS[j])
        ml = int(np.argmax(np.abs(acfs[j, 1:]))) + 1
        sig = sum(1 for lag in ASTRO_LAGS if abs(acfs[j, lag]) > ACF_THRESHOLD)
        fig3.suptitle(
            f"{LABELS[j]} ACF — biased circular autocorrelation, all 360 lags\n"
            f"max|ACF| = {max(abs(acfs[j, 1:])):.4f} @ lag {ml}°   ·   "
            f"95% CI = ±{ACF_THRESHOLD:.4f}   ·   "
            f"{sig}/{len(ASTRO_LAGS)} astrological lags significant   ·   "
            f"n = {n_planets} planets",
            fontsize=12, y=1.02
        )
        plt.tight_layout()
        out_sp = OUT_DIR / f"acf_{LABELS[j].lower()}.png"
        plt.savefig(out_sp, dpi=140, bbox_inches="tight")
        plt.close()
    print(f"  Per-SP plots → acf_<sp>.png  (11 files)")

    # ── CSV ──────────────────────────────────────────────────────── #
    with open(OUT_DIR / "acf_all_sps.csv", "w", newline="") as f:
        w = csv_mod.writer(f)
        header = ["lag_deg"] + [f"ACF_{LABELS[j]}" for j in range(n_sp)] + \
                  ["in_astro_lag"] + [f"sig_{LABELS[j]}" for j in range(n_sp)]
        w.writerow(header)
        for k in range(N_WAVE):
            row = [k] + [f"{acfs[j, k]:+.6e}" for j in range(n_sp)] + \
                  ["yes" if k in ASTRO_LAGS else "no"] + \
                  ["yes" if abs(acfs[j, k]) > ACF_THRESHOLD else "no"
                   for j in range(n_sp)]
            w.writerow(row)
    print(f"  CSV  → {OUT_DIR / 'acf_all_sps.csv'}")

    # ── Console table: ACF at astrological lags per SP ─────────── #
    print(f"\n══════ ACF at astrological lags, per SP (95% CI = ±{ACF_THRESHOLD:.4f}) ══════")
    hdr = f"  {'lag':<6} {'aspect':<10}  " + "  ".join(f"{LABELS[j][:7]:>8}" for j in range(n_sp))
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for lag in ASTRO_LAGS:
        row = f"  {lag:<6} {ASTRO_NAMES[lag]:<10}  "
        for j in range(n_sp):
            v = acfs[j, lag]
            marker = "*" if abs(v) > ACF_THRESHOLD else " "
            row += f"{v:>+7.3f}{marker}"
        print(row)
    print("  " + "-" * (len(hdr) - 2))
    print(f"  (*) = exceeds 95% white-noise CI")


if __name__ == "__main__":
    main()
