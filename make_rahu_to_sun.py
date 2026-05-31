"""
make_rahu_to_sun.py
    — Generate Rahu (mean north lunar node) Sun-relative longitude series.

  Produces a 290-element NumPy array of (Rahu_longitude − Sun_longitude) mod 360
  in degrees, for noon UTC on each of the 290 consecutive days from
  February 15, 2022 through December 1, 2022. This matches the encoding
  used by every other row of `corrected_sp_series.npy`.

  The mean lunar node longitude formula is from Meeus, *Astronomical
  Algorithms*, 2nd ed., Chapter 47:

      Ω = 125.0445550 − 1934.1361849·T + 0.0020762·T² + T³/467410 − T⁴/60616000

  where T is Julian centuries from J2000.0. The Sun's geometric mean
  longitude is from Meeus, Chapter 25.

  This script is provided for reproducibility. It also augments
  `corrected_sp_series.npy` to include Rahu as the 12th row when run.

  Run it directly:
      python make_rahu_to_sun.py
"""

from pathlib import Path
from datetime import datetime, timedelta
import numpy as np


OUT_DIR = Path(__file__).parent
SP_NPY  = OUT_DIR / "corrected_sp_series.npy"


def julian_date(dt: datetime) -> float:
    """Julian Date from a Python datetime treated as UTC."""
    y, m, d = dt.year, dt.month, dt.day
    h = dt.hour + dt.minute / 60.0 + dt.second / 3600.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5
    return jd + h / 24.0


def rahu_longitude(jd: float) -> float:
    """Mean longitude of the Moon's ascending node (Meeus Ch. 47)."""
    T = (jd - 2451545.0) / 36525.0
    omega = (
        125.0445550
        - 1934.1361849 * T
        + 0.0020762    * T**2
        + (T**3) / 467410.0
        - (T**4) / 60616000.0
    )
    return omega % 360.0


def sun_longitude(jd: float) -> float:
    """Sun's geometric mean longitude (Meeus Ch. 25)."""
    T = (jd - 2451545.0) / 36525.0
    L0 = (280.46646 + 36000.76983 * T + 0.0003032 * T**2) % 360.0
    M  = (357.52911 + 35999.05029 * T - 0.0001537 * T**2) % 360.0
    Mr = np.deg2rad(M)
    C = (
        (1.914602 - 0.004817 * T - 0.000014 * T**2) * np.sin(Mr)
        + (0.019993 - 0.000101 * T) * np.sin(2 * Mr)
        + 0.000289 * np.sin(3 * Mr)
    )
    return (L0 + C) % 360.0


def main():
    # 290 days, noon UTC, starting Feb 15 2022 — same grid as the paper
    dates = [datetime(2022, 2, 15, 12, 0) + timedelta(days=i) for i in range(290)]
    jds = np.array([julian_date(d) for d in dates])

    rahu_abs = np.array([rahu_longitude(jd) for jd in jds])
    sun_abs  = np.array([sun_longitude(jd)  for jd in jds])
    rahu_to_sun = (rahu_abs - sun_abs) % 360.0

    print(f"Date range: {dates[0]:%Y-%m-%d} to {dates[-1]:%Y-%m-%d}  ({len(dates)} days)")
    print(f"Rahu absolute  day 0: {rahu_abs[0]:.3f}°   day 289: {rahu_abs[-1]:.3f}°")
    print(f"Total Rahu motion:    {rahu_abs[-1] - rahu_abs[0]:+.3f}°   (expect ~ -15.3° = -19.3°/yr × 290/365)")
    print(f"Rahu-to-Sun    day 0: {rahu_to_sun[0]:.3f}°   day 289: {rahu_to_sun[-1]:.3f}°")

    # Save the standalone array
    out = OUT_DIR / "rahu_to_sun.npy"
    np.save(out, rahu_to_sun)
    print(f"\nSaved -> {out}")

    # Augment corrected_sp_series.npy if it has 11 rows
    if SP_NPY.exists():
        sp = np.load(SP_NPY)
        if sp.shape == (11, 290):
            sp_new = np.vstack([sp, rahu_to_sun[None, :]])
            np.save(SP_NPY, sp_new)
            print(f"Augmented {SP_NPY.name}: shape {sp.shape} -> {sp_new.shape}")
        elif sp.shape == (12, 290):
            # Already augmented; verify row 11 matches what we just computed
            mismatch = np.max(np.abs((sp[11] - rahu_to_sun + 180) % 360 - 180))
            print(f"{SP_NPY.name} already has 12 rows; row 11 max deviation from "
                  f"computed Rahu: {mismatch:.6f}°")
        else:
            print(f"Unexpected SP series shape {sp.shape}; not modifying")


if __name__ == "__main__":
    main()
