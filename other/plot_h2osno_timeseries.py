#!/usr/bin/env python3
"""Plot concatenated H2OSNO-derived time series from processed NetCDF output."""

import argparse
import os
from datetime import date, datetime

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset, num2date


def _read_1d_var(ds: Dataset, var_name: str) -> np.ndarray:
    if var_name not in ds.variables:
        raise KeyError(f"Variable '{var_name}' not found in {ds.filepath()}")
    arr = np.asarray(ds.variables[var_name][:]).squeeze()
    if arr.ndim != 1:
        raise ValueError(f"Variable '{var_name}' is not 1D after squeeze; shape={arr.shape}")
    return arr


def _is_leap_year(year: int, calendar: str) -> bool:
    if calendar in {"noleap", "365_day"}:
        return False
    if calendar in {"all_leap", "366_day"}:
        return True
    if calendar == "360_day":
        return False
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def _days_in_month(year: int, month: int, calendar: str) -> int:
    if calendar == "360_day":
        return 30

    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if month == 2 and _is_leap_year(year, calendar):
        return 29
    return month_lengths[month - 1]


def _dates_to_fractional_year(dates, calendar: str) -> np.ndarray:
    values = []
    for item in dates:
        year = item.year
        month = item.month
        day = item.day
        hour = getattr(item, "hour", 0)
        minute = getattr(item, "minute", 0)
        second = getattr(item, "second", 0)

        month_lengths = [_days_in_month(year, idx, calendar) for idx in range(1, 13)]
        days_before_month = sum(month_lengths[: month - 1])
        days_this_month = month_lengths[month - 1]
        day_fraction = ((day - 1) + hour / 24.0 + minute / 1440.0 + second / 86400.0) / days_this_month
        days_in_year = float(sum(month_lengths))
        year_fraction = (days_before_month + day_fraction * days_this_month) / days_in_year
        values.append(year + year_fraction)

    return np.asarray(values)


def _get_time_axis(ds: Dataset, n_time: int):
    if "time" not in ds.variables:
        return np.arange(n_time), "Time index"

    time_var = ds.variables["time"]
    time_vals = np.asarray(time_var[:]).squeeze()
    if time_vals.ndim != 1:
        return np.arange(n_time), "Time index"

    units = getattr(time_var, "units", None)
    cal = getattr(time_var, "calendar", "standard")
    if units is None:
        return np.arange(n_time), "Time index"

    try:
        dates = num2date(time_vals, units=units, calendar=cal)
        dates = np.atleast_1d(dates)
        if len(dates) == 0:
            return np.arange(n_time), "Time index"

        first = dates[0]
        if isinstance(first, (datetime, date)):
            return _dates_to_fractional_year(dates, cal), "Year"

        if first.__class__.__module__.startswith("cftime"):
            return _dates_to_fractional_year(dates, cal), "Year"

        return np.arange(n_time), "Time index"
    except Exception:
        return np.arange(n_time), "Time index"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot sumMaskedH2OSNOkm, mmSeaLevelEquiv, and mmSeaLevelEquivChange "
            "from a processed NetCDF file."
        )
    )
    parser.add_argument(
        "input_nc",
        help="Input NetCDF file created by sumH2OSNO.sh",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output image path (default: <input_stem>_timeseries.png)",
    )
    args = parser.parse_args()

    input_nc = args.input_nc
    if not os.path.isfile(input_nc):
        raise FileNotFoundError(f"Input file not found: {input_nc}")

    if args.output is None:
        stem, _ = os.path.splitext(input_nc)
        output_png = f"{stem}_timeseries.png"
    else:
        output_png = args.output

    with Dataset(input_nc, "r") as ds:
        v1 = _read_1d_var(ds, "sumMaskedH2OSNOkm")
        v2 = _read_1d_var(ds, "mmSeaLevelEquiv")
        v3 = _read_1d_var(ds, "mmSeaLevelEquivChange")

        n_time = len(v1)
        if len(v2) != n_time or len(v3) != n_time:
            raise ValueError("Expected all plotted variables to have the same length")

        x, xlabel = _get_time_axis(ds, n_time)

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True, constrained_layout=True)

    axes[0].plot(x, v1, color="tab:blue", lw=1.8)
    axes[0].set_ylabel("Gt")
    axes[0].set_title("liquid-water content of snowpack")
    axes[0].grid(alpha=0.3)

    axes[1].plot(x, v2, color="tab:green", lw=1.8)
    axes[1].set_ylabel("mm")
    axes[1].set_title("liquid-water content of snowpack (approx. sea-level equiv.)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(x, v3, color="tab:red", lw=1.8)
    axes[2].set_ylabel("mm")
    axes[2].set_title("liquid-water content of snowpack (approx. sea-level equiv. change)")
    axes[2].set_xlabel(xlabel)
    axes[2].grid(alpha=0.3)

    fig.suptitle("Greenland Snowpack-Derived Time Series", fontsize=14)
    fig.savefig(output_png, dpi=150)
    print(f"Wrote plot: {output_png}")
    plt.show()


if __name__ == "__main__":
    main()
