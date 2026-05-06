"""
plot_timeseries_ice.py
======================
Plot daily sea-ice metrics within a latitude band, comparing a snowcapping-fix
E3SM run against a control run.

Metrics plotted (one panel each):
  - Area-weighted mean sea-ice concentration within ±LATBAND_DEGREES
  - Maximum sea-ice concentration within ±LATBAND_DEGREES
  - Total sea-ice volume within ±LATBAND_DEGREES (km³)

Data source: *.mpassi.hist.am.timeSeriesStatsDaily.*.nc

Usage
-----
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_login.sh
    python plot_timeseries_ice.py

Edit config.py to change run directories, labels, or lat-band settings.
Output is saved to config.OUTPUT_DIR/ice_timeseries.png.
"""

import os
import sys

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Ensure config.py and utils.py in the same directory are importable.
_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)

import config
import utils


# ---------------------------------------------------------------------------
# Candidate variable names (tried in order; first match wins)
# ---------------------------------------------------------------------------

_ICE_AREA_CANDS = [
    'timeDaily_avg_iceAreaCell',
    'timeDaily_avg_iceArea',
    'iceAreaCell',
]
_ICE_VOL_CANDS = [
    'timeDaily_avg_iceVolumeCell',
    'timeDaily_avg_iceVolume',
    'iceVolumeCell',
]


def _find_var(ds, candidates, label):
    """Return the first candidate name found in ds, or None with a warning."""
    for name in candidates:
        if name in ds:
            return name
    avail = list(ds.data_vars)
    print(f"  Warning: no '{label}' variable found. "
          f"Tried: {candidates}. "
          f"Available (first 30): {avail[:30]}")
    return None


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_mpassi_daily(run_dir, latband_deg, mesh_lat, mesh_area):
    """
    Load mpassi timeSeriesStatsDaily files from run_dir and compute daily
    sea-ice metrics restricted to cells with |lat| ≤ latband_deg.

    Parameters
    ----------
    run_dir     : path to the simulation run directory
    latband_deg : latitude half-width for the analysis band
    mesh_lat    : 1-D array of cell latitudes in degrees (from restart file)
    mesh_area   : 1-D array of cell areas in m²       (from restart file)

    Returns a pd.DataFrame indexed by datetime with columns:
        mean_concentration  — area-weighted mean ice area fraction
        max_concentration   — maximum ice area fraction
        total_volume_km3    — total ice volume in km³

    If there is no ice within the lat band, metrics are zero (not NaN).
    Returns None if files are not found or cannot be opened.
    Works for multi-year runs (all matching files are loaded one at a time).

    Memory strategy: files are processed one at a time and only the in-band
    cell slice is kept in memory, avoiding OOM on login nodes.
    """
    files = utils.find_files(run_dir, '*.mpassi.hist.am.timeSeriesStatsDaily.*.nc')
    if not files:
        return None

    # --- Build latitude mask using mesh coords from restart ---------------
    in_band   = np.abs(mesh_lat) <= latband_deg
    n_band    = in_band.sum()
    area_band = mesh_area[in_band]
    print(f"    Cells in |lat| ≤ {latband_deg}°: {n_band} / {len(mesh_lat)}")

    # Probe first file to discover variable names
    try:
        ds0 = xr.open_dataset(files[0], decode_times=False, mask_and_scale=True)
    except Exception as exc:
        print(f"    Error opening {os.path.basename(files[0])}: {exc}")
        return None
    ia_var = _find_var(ds0, _ICE_AREA_CANDS, 'ice area fraction')
    iv_var = _find_var(ds0, _ICE_VOL_CANDS,  'ice volume')
    ds0.close()

    if ia_var is None and iv_var is None:
        print("    No ice area or volume variables found; nothing to plot.")
        return None

    # --- Process one file at a time to keep memory low -------------------
    all_times  = []
    mean_conc  = []
    max_conc   = []
    total_vol  = []

    print(f"    Processing {len(files)} file(s) one at a time…")
    for fpath in files:
        try:
            ds = xr.open_dataset(fpath, decode_times=False, mask_and_scale=True)
        except Exception as exc:
            print(f"    Warning: cannot open {os.path.basename(fpath)}: {exc}")
            continue

        n_t = ds.sizes.get('Time', 0)

        # Time axis
        t = utils.parse_mpas_xtime(ds)
        if t is not None and len(t) == n_t:
            all_times.extend(t)
        else:
            # Fall back: integer placeholders (replaced later with a range index)
            all_times.extend([None] * n_t)

        # Ice area
        if ia_var is not None and ia_var in ds:
            ia = ds[ia_var].values[:, in_band].astype(float)  # (Time, nBand)
            ia = np.where(np.abs(ia) < 1e30, ia, np.nan)

            wt     = np.where(np.isfinite(ia), area_band[np.newaxis, :], 0.0)
            wt_sum = wt.sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                mc = np.where(wt_sum > 0,
                              np.nansum(ia * area_band[np.newaxis, :], axis=1) / wt_sum,
                              0.0)
            mean_conc.extend(mc.tolist())

            with np.errstate(all='ignore'):
                mx = np.nanmax(ia, axis=1)
            max_conc.extend(np.where(np.isfinite(mx), mx, 0.0).tolist())
        else:
            mean_conc.extend([np.nan] * n_t)
            max_conc.extend([np.nan] * n_t)

        # Ice volume
        if iv_var is not None and iv_var in ds:
            iv = ds[iv_var].values[:, in_band].astype(float)
            iv = np.where(np.abs(iv) < 1e30, iv, np.nan)
            total_vol.extend(
                (np.nansum(iv * area_band[np.newaxis, :], axis=1) / 1e9).tolist()
            )
        else:
            total_vol.extend([np.nan] * n_t)

        ds.close()

    # --- Assemble DataFrame -----------------------------------------------
    results = {}
    if ia_var is not None:
        results['mean_concentration'] = mean_conc
        results['max_concentration']  = max_conc
    if iv_var is not None:
        results['total_volume_km3'] = total_vol

    df = pd.DataFrame(results)

    # Build time index — use parsed datetimes if all were recovered
    if all(t is not None for t in all_times):
        df.index = pd.DatetimeIndex(all_times)
    else:
        df.index = pd.RangeIndex(len(df))
        print("    Warning: partial time axis; using integer index.")

    print(f"    Loaded {len(df)} daily records.")
    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

_PANEL_META = {
    'mean_concentration': {
        'title': 'Area-weighted Mean Sea-Ice Concentration',
        'ylabel': 'Ice fraction (0–1)',
    },
    'max_concentration': {
        'title': 'Maximum Sea-Ice Concentration',
        'ylabel': 'Ice fraction (0–1)',
    },
    'total_volume_km3': {
        'title': 'Total Sea-Ice Volume',
        'ylabel': 'Volume (km³)',
    },
}


def make_plot(fix_df, ctrl_df, latband_deg, output_dir):
    """
    Save ice_timeseries.png to output_dir.
    fix_df and ctrl_df may be None (the plot skips whichever is missing).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Determine which panels to draw
    all_cols = set()
    if fix_df  is not None: all_cols |= set(fix_df.columns)
    if ctrl_df is not None: all_cols |= set(ctrl_df.columns)
    cols = [c for c in _PANEL_META if c in all_cols]

    if not cols:
        print("  No plottable columns found for ice time series; skipping.")
        return

    fig, axes = plt.subplots(len(cols), 1,
                             figsize=(13, 3.2 * len(cols)),
                             sharex=True)
    if len(cols) == 1:
        axes = [axes]

    for ax, col in zip(axes, cols):
        meta = _PANEL_META[col]
        if fix_df is not None and col in fix_df.columns:
            ax.plot(fix_df.index, fix_df[col],
                    color=config.COLORS['fix'],
                    ls=config.LINESTYLES['fix'],
                    lw=config.LINEWIDTHS['fix'],
                    label=config.FIX_LABEL, zorder=3)
        if ctrl_df is not None and col in ctrl_df.columns:
            ax.plot(ctrl_df.index, ctrl_df[col],
                    color=config.COLORS['ctrl'],
                    ls=config.LINESTYLES['ctrl'],
                    lw=config.LINEWIDTHS['ctrl'],
                    label=config.CTRL_LABEL, zorder=2)
        ax.set_ylabel(meta['ylabel'], fontsize=9)
        ax.set_title(f"{meta['title']}  (|lat| ≤ {latband_deg:.0f}°)",
                     fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, linewidth=0.5)

    # Format x-axis as months when a DatetimeIndex is available
    ref_idx = (fix_df.index if fix_df is not None and len(fix_df) > 0
               else ctrl_df.index if ctrl_df is not None else None)
    if ref_idx is not None and isinstance(ref_idx, pd.DatetimeIndex):
        axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
        plt.setp(axes[-1].get_xticklabels(), rotation=30, ha='right')

    axes[-1].set_xlabel('Date')
    fig.suptitle('Sea-Ice Metrics Comparison within Latitude Band',
                 fontsize=12, y=1.01)
    fig.tight_layout()

    outpath = os.path.join(output_dir, 'ice_timeseries.png')
    fig.savefig(outpath, dpi=config.DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Sea-Ice Time Series")
    print("=" * 60)
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    latband = utils.resolve_latband(config.LATBAND_DEGREES, config.FIX_RUN_DIR)

    # Load mesh coordinates once from the restart file of each run.
    # The two simulations share the same MPAS mesh, so either restart works;
    # we load from each run independently in case they differ in future use.
    print(f"\nLoading mesh from fix run…")
    try:
        fix_lat, _fix_lon, fix_area = utils.load_mpassi_mesh(config.FIX_RUN_DIR)
    except RuntimeError as e:
        print(f"  Error: {e}")
        fix_lat = fix_area = None

    print(f"\nLoading mesh from ctrl run…")
    try:
        ctrl_lat, _ctrl_lon, ctrl_area = utils.load_mpassi_mesh(config.CTRL_RUN_DIR)
    except RuntimeError as e:
        print(f"  Error: {e}")
        ctrl_lat = ctrl_area = None

    print(f"\nLoading fix run data:  {config.FIX_RUN_DIR}")
    fix_df  = (load_mpassi_daily(config.FIX_RUN_DIR,  latband, fix_lat,  fix_area)
               if fix_lat is not None else None)

    print(f"\nLoading ctrl run data: {config.CTRL_RUN_DIR}")
    ctrl_df = (load_mpassi_daily(config.CTRL_RUN_DIR, latband, ctrl_lat, ctrl_area)
               if ctrl_lat is not None else None)

    print("\nGenerating plot…")
    make_plot(fix_df, ctrl_df, latband, config.OUTPUT_DIR)
    print("Done.\n")


if __name__ == '__main__':
    main()
