"""
effect_size_plot.py
    — Figure 6: per-reference-point effect-size dot plot.

  Reads acf_effect_size_results.csv (produced by acf_effect_size.py) and
  renders the standardised effect d for every reference point against both
  null models -- the compound null and the unnamed-MBA (named-vs-unnamed)
  null. Each reference point is a row: a filled navy marker shows d
  against the compound null, a hollow orange marker shows d against the
  unnamed-MBA null, and a connector between them makes the attenuation
  visible at a glance.

  Three one-sided significance thresholds (p < .05, .01, .001) are drawn
  as dashed verticals at d = 1.645, 2.326, 3.090. The Moon, whose mean
  d is negative under both nulls, is visually separated as a negative
  control. The omnibus row that combines all eleven points is drawn at
  the bottom with diamond markers.

  The RMS autocorrelation and Cohen-style magnitude label for each row
  are printed in a right-hand annotation column.

  Output: effect_size_dotplot.png. Run it directly:
      python effect_size_plot.py
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUT_DIR  = Path(__file__).parent
CSV_FILE = OUT_DIR / "acf_effect_size_results.csv"
OUT_PNG  = OUT_DIR / "effect_size_dotplot.png"

# Colours and styles
BG       = "#f4f0e6"     # warm cream
COMPOUND = "#1f3a93"     # navy
UNNAMED  = "#d99500"     # gold/orange
GRID     = "#9aa0a6"
MOON_RED = "#8a2b2b"

# One-sided z-thresholds and labels
THRESHOLDS = [(1.645, "p < .05"),
              (2.326, "p < .01"),
              (3.090, "p < .001")]


def main():
    # ---- read CSV ------------------------------------------------------ #
    rows = []
    with open(CSV_FILE) as f:
        for r in csv.DictReader(f):
            rows.append(r)

    summary = None
    sp_rows = []
    for r in rows:
        if r["special_point"].upper().startswith("ALL"):
            summary = r
        else:
            sp_rows.append(r)

    # Separate the Moon as a negative control if its d_compound is negative.
    moon_row = None
    main_rows = []
    for r in sp_rows:
        if r["special_point"] == "Moon" and float(r["d_compound"]) < 0:
            moon_row = r
        else:
            main_rows.append(r)

    # Sort the main rows by d_compound descending so the strongest sits at the
    # top of the figure.
    main_rows.sort(key=lambda r: float(r["d_compound"]), reverse=True)

    # Final row order, top to bottom in the figure:
    #   main_rows (sorted), then Moon (if separated), then omnibus.
    plot_rows = main_rows + ([moon_row] if moon_row is not None else [])
    n_main = len(plot_rows)
    n_total = n_main + (1 if summary is not None else 0)

    # y positions: top of figure = first row; omnibus at the very bottom
    # below a small visual gap.
    y_main = np.arange(n_main, 0, -1, dtype=float)        # n_main..1
    y_omni = 0.0
    y_gap = -0.5                                          # divider line y

    # ---- figure -------------------------------------------------------- #
    fig, ax = plt.subplots(figsize=(13, 7.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor(BG)

    # Significance threshold verticals
    for x, lab in THRESHOLDS:
        ax.axvline(x, color=GRID, lw=1.0, ls="--", alpha=0.7, zorder=1)
        ax.text(x, n_main + 1.05, lab, ha="center", va="bottom",
                fontsize=10, color="#555")
    ax.axvline(0, color="black", lw=1.0, alpha=0.9, zorder=1)

    # Connectors and markers
    for yi, r in zip(y_main, plot_rows):
        dc = float(r["d_compound"])
        dm = float(r["d_unnamed_mba"])
        ax.plot([dc, dm], [yi, yi], color=GRID, lw=1.2, alpha=0.7, zorder=2)
        ax.scatter([dc], [yi], s=130, marker="o", facecolor=COMPOUND,
                   edgecolor=COMPOUND, zorder=4)
        ax.scatter([dm], [yi], s=130, marker="o", facecolor="white",
                   edgecolor=UNNAMED, linewidth=2.0, zorder=4)

    # Omnibus row
    if summary is not None:
        dc_s = float(summary["d_compound"])
        dm_s = float(summary["d_unnamed_mba"])
        ax.plot([dc_s, dm_s], [y_omni, y_omni], color=GRID, lw=1.4, alpha=0.7,
                zorder=2)
        ax.scatter([dc_s], [y_omni], s=240, marker="D", facecolor=COMPOUND,
                   edgecolor=COMPOUND, zorder=4)
        ax.scatter([dm_s], [y_omni], s=240, marker="D", facecolor="white",
                   edgecolor=UNNAMED, linewidth=2.2, zorder=4)
        ax.axhline(y_gap, color=GRID, lw=0.6, alpha=0.6)

    # Y tick labels with Moon shown in italic red as a negative control
    yticks = list(y_main) + ([y_omni] if summary is not None else [])
    yticklabels = []
    label_colors = []
    label_styles = []
    for r in plot_rows:
        if r["special_point"] == "Moon":
            yticklabels.append("Moon  (neg. control)")
            label_colors.append(MOON_RED)
            label_styles.append("italic")
        else:
            yticklabels.append(r["special_point"])
            label_colors.append("black")
            label_styles.append("normal")
    if summary is not None:
        yticklabels.append("Omnibus (all 11)")
        label_colors.append("black")
        label_styles.append("normal")

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=11)
    # Style individual y-tick labels
    for tick, col, style in zip(ax.get_yticklabels(), label_colors, label_styles):
        tick.set_color(col)
        tick.set_fontstyle(style)
        if "Omnibus" in tick.get_text():
            tick.set_fontweight("bold")

    # X-axis limits
    all_d = [float(r["d_compound"])    for r in sp_rows] + \
            [float(r["d_unnamed_mba"]) for r in sp_rows]
    if summary is not None:
        all_d += [float(summary["d_compound"]), float(summary["d_unnamed_mba"])]
    xmin = min(all_d) - 0.7
    xmax = max(all_d) + 0.7
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(y_gap - 0.6, n_main + 1.6)

    # Right-hand RMS annotation column
    x_rms = xmax + (xmax - xmin) * 0.025
    ax.text(x_rms, n_main + 1.05, "RMS  r̄",
            fontsize=11, fontweight="bold", va="bottom")
    ax.text(x_rms, n_main + 0.55, "(abs. effect)",
            fontsize=9, style="italic", color="#666", va="bottom")
    for yi, r in zip(y_main, plot_rows):
        ax.text(x_rms, yi, f"{float(r['rms_acf']):.3f}",
                fontsize=10, va="center")
        ax.text(x_rms + (xmax - xmin) * 0.08, yi, r["magnitude"],
                fontsize=10, style="italic", color="#555", va="center")
    if summary is not None:
        ax.text(x_rms, y_omni, f"{float(summary['rms_acf']):.3f}",
                fontsize=11, fontweight="bold", va="center")
        ax.text(x_rms + (xmax - xmin) * 0.08, y_omni, summary["magnitude"],
                fontsize=10, style="italic", color="#555", va="center")

    # Legend (top-left inside the axes)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=11,
               markerfacecolor=COMPOUND, markeredgecolor=COMPOUND,
               label="Compound null"),
        Line2D([0], [0], marker="o", linestyle="", markersize=11,
               markerfacecolor="white", markeredgecolor=UNNAMED,
               markeredgewidth=2.0, label="Unnamed-asteroid null"),
    ]
    ax.legend(handles=handles, loc="upper left", framealpha=0.95,
              fontsize=10, edgecolor="#cccccc")

    # Titles
    ax.set_xlabel("Standardised effect  d   (null standard deviations above chance)",
                  fontsize=12, labelpad=8)
    fig.suptitle("Effect size of each reference point, tested against two null models",
                 fontsize=14, fontweight="bold", y=0.99, x=0.5)
    ax.set_title(
        "One-sided tests. Points right of a dashed line clear that "
        "significance level. Each line shows the drop from the compound "
        "null to the stricter unnamed-asteroid (named-vs-unnamed) null.",
        fontsize=10, style="italic", color="#444", pad=20)

    # Cosmetics
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors="#444")

    plt.subplots_adjust(left=0.12, right=0.82, top=0.86, bottom=0.10)
    plt.savefig(OUT_PNG, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved -> {OUT_PNG}")


if __name__ == "__main__":
    main()
