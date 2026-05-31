"""
heatmap_acf_normalized_by_joint.py
   — Heatmap of  ACF_lag_X / joint_acf_sq_obs  per SP, separating
     classical BENEFICS (Jupiter, Venus) from classical MALEFICS
     (Saturn, Mars).

  Why this normalisation
  ----------------------
  joint_acf_sq_obs = Σ ACF(ℓ)²  over the 11 astrological lags.  Dividing
  each ACF(ℓ) value by this total gives a per-SP normalised measure of
  *which lag dominates the SP's joint signal*.  Cells with large
  magnitude → that lag contributes a disproportionate fraction of the
  SP's overall astrological correlation budget.

  Note: division is by joint_acf_sq (a sum of squares), not by the
  L₂-norm √(joint_acf_sq); this is the literal "ACF / joint_acf_sq"
  the user requested.  A complementary view  ACF² / joint_acf_sq  is
  also produced — that is the variance-share interpretation and sums
  to 1 across the 11 lags for each SP.

  Outputs
  -------
  1. Three-row figure: BENEFICS heatmap, MALEFICS heatmap, Δ (ben − mal)
     row — all using ACF / joint_acf_sq.
  2. Same layout using the variance-share ACF² / joint_acf_sq view.
  3. Full 11-SP grouped-by-category heatmap (ACF / joint_acf_sq).
  4. CSV of normalised values.
"""

import csv as csv_mod
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT_DIR  = Path(__file__).parent
CSV_FILE = OUT_DIR / "acf_significance_results.csv"

N_WAVE = 360
ACF_THRESHOLD = 1.96 / np.sqrt(N_WAVE)

ASTRO_LAGS = [24, 30, 36, 40, 45, 51, 60, 72, 90, 120, 180]
ASTRO_NAMES = {
    24:  "quindecile",   30:  "semi-sextile", 36:  "decile",
    40:  "novile",       45:  "semi-square",  51:  "septile",
    60:  "sextile",      72:  "quintile",     90:  "square",
    120: "trine",        180: "opposition",
}

GROUP_ORDER = [
    ("Classical BENEFICS",   ["Jupiter", "Venus"],                "#1b7a3a"),
    ("Classical MALEFICS",   ["Saturn",  "Mars"],                 "#a01717"),
    ("Luminaries",           ["Sun",     "Moon"],                 "#d68910"),
    ("Neutral",              ["Mercury"],                          "#7a5cb1"),
    ("Reference",            ["Mesarthim"],                        "#666666"),
    ("Modern (post-class.)", ["Uranus",  "Neptune", "Pluto"],     "#2c5896"),
    ("Lunar Node (Vedic)",   ["Rahu"],                             "#6a1b9a"),
]


def load_acf_results(path):
    rows = {}
    with open(path) as f:
        reader = csv_mod.DictReader(f)
        for row in reader:
            if not row.get("special_point"):
                continue
            rows[row["special_point"]] = {
                "K":      int(row["K_acf_obs"]),
                "JOINT":  float(row["joint_acf_sq_obs"]),
                "LBQ":    float(row["LB_Q_obs"]),
                "acf":    np.array([float(row[f"ACF_lag_{lag}"]) for lag in ASTRO_LAGS]),
            }
    return rows


def annotate_cell(ax, ri, ci, value, font_size=9, threshold=0.20,
                   fmt="{:+.2f}"):
    color = "white" if abs(value) > threshold else "black"
    ax.text(ci, ri, fmt.format(value), ha="center", va="center",
             fontsize=font_size, color=color)


def main():
    rows = load_acf_results(CSV_FILE)
    print(f"Loaded ACF significance results for {len(rows)} SPs.")
    print(f"\nPer-SP joint_acf_sq_obs (denominator):")
    for sp, d in rows.items():
        print(f"  {sp:<12} joint_acf_sq = {d['JOINT']:.5f}   "
              f"K = {d['K']:<2d}  LB_Q = {d['LBQ']:.2f}")

    col_labels = [f"{lag}°\n{ASTRO_NAMES[lag]}" for lag in ASTRO_LAGS]

    # Compute normalised matrices
    benefics = ["Jupiter", "Venus"]
    malefics = ["Saturn",  "Mars"]

    def norm_acf(sp_list):
        out = []
        for sp in sp_list:
            out.append(rows[sp]["acf"] / rows[sp]["JOINT"])
        return np.array(out)

    def norm_acf_sq(sp_list):
        out = []
        for sp in sp_list:
            out.append(rows[sp]["acf"]**2 / rows[sp]["JOINT"])
        return np.array(out)

    ben_norm   = norm_acf(benefics)        # ACF / joint
    mal_norm   = norm_acf(malefics)
    ben_share  = norm_acf_sq(benefics)     # ACF² / joint
    mal_share  = norm_acf_sq(malefics)

    diff_norm  = (ben_norm.mean(axis=0)  - mal_norm.mean(axis=0))[None, :]
    diff_share = (ben_share.mean(axis=0) - mal_share.mean(axis=0))[None, :]

    # ── FIGURE 1: ACF / joint_acf_sq — benefics vs malefics ────── #
    vlim = float(max(np.max(np.abs(ben_norm)), np.max(np.abs(mal_norm)))) * 1.05
    fig, axes = plt.subplots(3, 1, figsize=(15, 9),
                              gridspec_kw=dict(height_ratios=[2, 2, 1.2], hspace=0.55))

    for ax, mat, names, title in [
        (axes[0], ben_norm, benefics,
            "Classical BENEFICS  —  ACF / joint_acf_sq_obs"),
        (axes[1], mal_norm, malefics,
            "Classical MALEFICS  —  ACF / joint_acf_sq_obs"),
    ]:
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vlim, vmax=+vlim, aspect="auto")
        for ri in range(mat.shape[0]):
            for ci in range(mat.shape[1]):
                v = mat[ri, ci]
                color = "white" if abs(v) > vlim * 0.55 else "black"
                ax.text(ci, ri, f"{v:+.2f}", ha="center", va="center",
                         fontsize=10, color=color)
                # Outline cells where the underlying ACF exceeds 95% CI
                sp_idx = (benefics + malefics).index(names[ri]) if names[ri] in (benefics + malefics) else None
                raw_acf = rows[names[ri]]["acf"][ci]
                if abs(raw_acf) > ACF_THRESHOLD:
                    ax.add_patch(Rectangle((ci - 0.5, ri - 0.5), 1, 1,
                                            linewidth=1.6, edgecolor="black",
                                            facecolor="none"))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels([f"{n}\n(jt={rows[n]['JOINT']:.3f})" for n in names], fontsize=10)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
        ax.set_title(title, fontsize=12, pad=8)
        fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="ACF / joint_acf_sq")

    # Δ row
    dmax = float(np.max(np.abs(diff_norm))) * 1.1
    im_d = axes[2].imshow(diff_norm, cmap="PiYG", vmin=-dmax, vmax=+dmax,
                           aspect="auto")
    for ci in range(diff_norm.shape[1]):
        v = diff_norm[0, ci]
        color = "white" if abs(v) > dmax * 0.55 else "black"
        axes[2].text(ci, 0, f"{v:+.2f}", ha="center", va="center",
                      fontsize=11, color=color, fontweight="bold")
    axes[2].set_xticks(np.arange(diff_norm.shape[1]))
    axes[2].set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    axes[2].set_yticks([0])
    axes[2].set_yticklabels(["⟨BENEFICS⟩ − ⟨MALEFICS⟩"], fontsize=10,
                              fontweight="bold")
    axes[2].set_title(
        "Δ = mean(BENEFICS) − mean(MALEFICS)  of  ACF / joint_acf_sq\n"
        "(positive ⇒ benefics carry more of their joint signal at that aspect)",
        fontsize=11, pad=8
    )
    fig.colorbar(im_d, ax=axes[2], fraction=0.025, pad=0.01, label="Δ")

    fig.suptitle(
        "Per-SP-normalised ACF heatmap:  ACF(ℓ) / joint_acf_sq_obs\n"
        f"Black outline = |ACF(ℓ)| > 95% CI band (±{ACF_THRESHOLD:.4f}).  "
        f"Each row's normalisation = sum of squared ACF over all 11 astrological lags.",
        fontsize=12, y=1.005
    )
    plt.tight_layout()
    out1 = OUT_DIR / "heatmap_acf_norm_by_joint_benefics_vs_malefics.png"
    plt.savefig(out1, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot → {out1}")

    # ── FIGURE 2: variance-share view (ACF² / joint) ─────────────── #
    fig2, axes2 = plt.subplots(3, 1, figsize=(15, 9),
                                gridspec_kw=dict(height_ratios=[2, 2, 1.2], hspace=0.55))
    share_max = float(max(np.max(ben_share), np.max(mal_share)))
    for ax, mat, names, title in [
        (axes2[0], ben_share, benefics,
            "Classical BENEFICS  —  ACF² / joint_acf_sq  (variance share, sums to 1 per row)"),
        (axes2[1], mal_share, malefics,
            "Classical MALEFICS  —  ACF² / joint_acf_sq  (variance share, sums to 1 per row)"),
    ]:
        im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=share_max, aspect="auto")
        for ri in range(mat.shape[0]):
            for ci in range(mat.shape[1]):
                v = mat[ri, ci]
                color = "white" if v < share_max * 0.45 else "black"
                ax.text(ci, ri, f"{v*100:.1f}%", ha="center", va="center",
                         fontsize=10, color=color)
                raw_acf = rows[names[ri]]["acf"][ci]
                if abs(raw_acf) > ACF_THRESHOLD:
                    ax.add_patch(Rectangle((ci - 0.5, ri - 0.5), 1, 1,
                                            linewidth=1.6, edgecolor="orange",
                                            facecolor="none"))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels([f"{n}\n(Σ rows = 100%)" for n in names], fontsize=10)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
        ax.set_title(title, fontsize=12, pad=8)
        fig2.colorbar(im, ax=ax, fraction=0.025, pad=0.01,
                       label="ACF² / joint  (fraction)")
    # Δ row of variance share
    dmax2 = float(np.max(np.abs(diff_share))) * 1.1
    im_d2 = axes2[2].imshow(diff_share, cmap="PiYG", vmin=-dmax2, vmax=+dmax2,
                              aspect="auto")
    for ci in range(diff_share.shape[1]):
        v = diff_share[0, ci]
        color = "white" if abs(v) > dmax2 * 0.55 else "black"
        axes2[2].text(ci, 0, f"{v*100:+.1f}%", ha="center", va="center",
                       fontsize=11, color=color, fontweight="bold")
    axes2[2].set_xticks(np.arange(diff_share.shape[1]))
    axes2[2].set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    axes2[2].set_yticks([0])
    axes2[2].set_yticklabels(["⟨BENEFICS⟩ − ⟨MALEFICS⟩"], fontsize=10,
                              fontweight="bold")
    axes2[2].set_title(
        "Δ variance share  =  mean(BENEFICS) − mean(MALEFICS)\n"
        "(positive ⇒ benefics allocate a larger fraction of their joint signal to that aspect)",
        fontsize=11, pad=8
    )
    fig2.colorbar(im_d2, ax=axes2[2], fraction=0.025, pad=0.01, label="Δ share")
    fig2.suptitle(
        "Variance-share heatmap:  ACF(ℓ)² / joint_acf_sq_obs   "
        "(each row sums to 100% across the 11 astrological lags)\n"
        "Orange outline = |ACF(ℓ)| > 95% CI band.",
        fontsize=12, y=1.005
    )
    plt.tight_layout()
    out2 = OUT_DIR / "heatmap_acf_variance_share_benefics_vs_malefics.png"
    plt.savefig(out2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot → {out2}")

    # ── FIGURE 3: full 11-SP grouped heatmap, ACF/joint ────────── #
    ordered_sps = []
    group_separators = []
    group_color_bar = []
    cursor = 0
    for group_name, group_sps, group_color in GROUP_ORDER:
        present = [sp for sp in group_sps if sp in rows]
        if not present: continue
        if cursor > 0:
            group_separators.append(cursor)
        start = cursor
        for sp in present:
            ordered_sps.append(sp)
            cursor += 1
        end = cursor - 1
        group_color_bar.append((start, end, group_color, group_name))

    mat_all = np.array([rows[sp]["acf"] / rows[sp]["JOINT"] for sp in ordered_sps])
    sig_mask = np.array([np.abs(rows[sp]["acf"]) > ACF_THRESHOLD for sp in ordered_sps])

    vlim_all = float(np.max(np.abs(mat_all))) * 1.05
    fig3 = plt.figure(figsize=(17, 10))
    gs = fig3.add_gridspec(1, 2, width_ratios=[14, 1], wspace=0.05)
    ax_m = fig3.add_subplot(gs[0, 0])
    ax_c = fig3.add_subplot(gs[0, 1])
    im = ax_m.imshow(mat_all, cmap="RdBu_r", vmin=-vlim_all, vmax=+vlim_all,
                       aspect="auto")
    for ri in range(mat_all.shape[0]):
        for ci in range(mat_all.shape[1]):
            v = mat_all[ri, ci]
            color = "white" if abs(v) > vlim_all * 0.55 else "black"
            ax_m.text(ci, ri, f"{v:+.2f}", ha="center", va="center",
                       fontsize=9, color=color)
            if sig_mask[ri, ci]:
                ax_m.add_patch(Rectangle((ci - 0.5, ri - 0.5), 1, 1,
                                          linewidth=1.4, edgecolor="black",
                                          facecolor="none"))
    ax_m.set_xticks(np.arange(mat_all.shape[1]))
    ax_m.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax_m.set_yticks(np.arange(mat_all.shape[0]))
    ax_m.set_yticklabels(
        [f"{sp}\n(jt={rows[sp]['JOINT']:.3f})" for sp in ordered_sps], fontsize=10
    )
    for sep in group_separators:
        ax_m.axhline(sep - 0.5, color="black", lw=1.5)
    for (start, end, color, label) in group_color_bar:
        ax_m.add_patch(Rectangle((-2.2, start - 0.5), 0.6, (end - start + 1),
                                  facecolor=color, edgecolor="black",
                                  linewidth=0.7, clip_on=False))
        ax_m.text(-2.4, (start + end) / 2, label, ha="right",
                  va="center", fontsize=9, fontweight="bold", color=color)
    ax_m.set_xlabel("Astrological aspect", fontsize=11)
    ax_m.set_title(
        "All 11 SPs (grouped by classical category)  —  ACF / joint_acf_sq_obs\n"
        f"black outline = |ACF| > 95% CI (±{ACF_THRESHOLD:.4f})",
        fontsize=12, pad=10
    )
    fig3.colorbar(im, cax=ax_c, label="ACF / joint")
    plt.tight_layout()
    out3 = OUT_DIR / "heatmap_acf_norm_by_joint_all_grouped.png"
    plt.savefig(out3, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Plot → {out3}")

    # ── CSV ─────────────────────────────────────────────────────── #
    with open(OUT_DIR / "acf_norm_by_joint.csv", "w", newline="") as f:
        w = csv_mod.writer(f)
        header = ["special_point", "category", "joint_acf_sq_obs"] + \
                 [f"acf_norm_lag_{lag}" for lag in ASTRO_LAGS] + \
                 [f"variance_share_lag_{lag}" for lag in ASTRO_LAGS]
        w.writerow(header)
        cat_lookup = {sp: grp for grp, sps, _ in GROUP_ORDER for sp in sps}
        for sp in ordered_sps:
            d = rows[sp]
            cat = cat_lookup.get(sp, "")
            row = [sp, cat, f"{d['JOINT']:.6f}"]
            for i, lag in enumerate(ASTRO_LAGS):
                row.append(f"{d['acf'][i] / d['JOINT']:+.6f}")
            for i, lag in enumerate(ASTRO_LAGS):
                row.append(f"{d['acf'][i]**2 / d['JOINT']:.6f}")
            w.writerow(row)
    print(f"  CSV  → {OUT_DIR / 'acf_norm_by_joint.csv'}")

    # ── Console tables ─────────────────────────────────────────── #
    print(f"\n══════ ACF / joint_acf_sq_obs ──  Benefics & Malefics ══════")
    print(f"  {'lag':<6} {'aspect':<14} {'Jupiter':>8} {'Venus':>8} {'⟨ben⟩':>8} | "
          f"{'Saturn':>8} {'Mars':>8} {'⟨mal⟩':>8} | {'Δ ben−mal':>10}")
    print("  " + "-" * 95)
    for i, lag in enumerate(ASTRO_LAGS):
        j_n  = rows["Jupiter"]["acf"][i] / rows["Jupiter"]["JOINT"]
        v_n  = rows["Venus"]["acf"][i]   / rows["Venus"]["JOINT"]
        s_n  = rows["Saturn"]["acf"][i]  / rows["Saturn"]["JOINT"]
        m_n  = rows["Mars"]["acf"][i]    / rows["Mars"]["JOINT"]
        bm   = (j_n + v_n) / 2; mm = (s_n + m_n) / 2
        print(f"  {lag:<6} {ASTRO_NAMES[lag]:<14} {j_n:>+8.3f} {v_n:>+8.3f} "
              f"{bm:>+8.3f} | {s_n:>+8.3f} {m_n:>+8.3f} {mm:>+8.3f} | "
              f"{(bm - mm):>+10.3f}")

    print(f"\n══════ Variance share  ACF² / joint  (each row sums to 100%) ══════")
    print(f"  {'lag':<6} {'aspect':<14} {'Jupiter':>8} {'Venus':>8} {'⟨ben⟩':>8} | "
          f"{'Saturn':>8} {'Mars':>8} {'⟨mal⟩':>8} | {'Δ ben−mal':>10}")
    print("  " + "-" * 95)
    for i, lag in enumerate(ASTRO_LAGS):
        j_s  = rows["Jupiter"]["acf"][i]**2 / rows["Jupiter"]["JOINT"]
        v_s  = rows["Venus"]["acf"][i]**2   / rows["Venus"]["JOINT"]
        s_s  = rows["Saturn"]["acf"][i]**2  / rows["Saturn"]["JOINT"]
        m_s  = rows["Mars"]["acf"][i]**2    / rows["Mars"]["JOINT"]
        bm   = (j_s + v_s) / 2; mm = (s_s + m_s) / 2
        print(f"  {lag:<6} {ASTRO_NAMES[lag]:<14} {j_s*100:>7.2f}% {v_s*100:>7.2f}% "
              f"{bm*100:>7.2f}% | {s_s*100:>7.2f}% {m_s*100:>7.2f}% {mm*100:>7.2f}% | "
              f"{(bm - mm)*100:>+9.2f}%")


if __name__ == "__main__":
    main()
