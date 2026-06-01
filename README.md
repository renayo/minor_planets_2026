# Minor planet orbits and astronyms in the news

This repository contains the full analysis pipeline, data, and figures for
Oshop & Coops, *In daily news, minor planet names show up in accordance with harmonic patterning* (2026).

The study tracks the proper names of 1,121 officially designated minor
planets across Google News for 290 consecutive days in 2022, then asks
whether the angular separations between named minor planets and twelve
reference points show statistical structure at the harmonic angles long
treated as salient in astrological tradition.

## Reference points

Twelve reference points are evaluated:

1. **Mesarthim** — fixed star at γ Arietis, also serves as a proxy for the
   First Point of Aries (vernal equinox; the values used are absolute
   tropical longitudes near 0° Aries).
2. **Sun** — solar tropical longitude.
3. **Moon** — lunar tropical longitude.
4. **Mercury** through **Pluto** — geocentric tropical longitudes of the
   seven classical and modern planets (5 through 10).
5. **Rahu**: the mean north lunar node. A
   mathematical point central to Vedic astrology, computed from the
   standard Meeus polynomial for the mean longitude of the ascending node
   (Astronomical Algorithms, Chapter 47).

## Quick-start reproduction

See `HOW TO REPLICATE.md` for a step-by-step guide. The main analysis
scripts can be run in any order after the data and `corrected_sp_series.npy`
(12 × 290 array) are in place.

```bash
python3 acf_significance.py        # per-SP and omnibus ACF significance
python3 "all v test.py"            # V-test at every direction, two nulls
python3 acf_effect_size.py         # standardised effect d, per-SP and omnibus
python3 acf_plots_all_sps.py       # Figure 5 ACF plots (linear + polar)
python3 effet_size_plot.py         # Figure 6 lollipop effect-size chart
python3 all_sps_vs_gauquelin_extracted.py   # Figures 7-8 Gauquelin overlay
python3 heatmap_acf_normalized_by_joint.py  # benefic/malefic heatmaps
```

The Monte Carlo iteration count was set to **5,000**.

## Paper

Under submission.
