"""
plot_timeseries_runoff.py
=========================
Plot monthly time series of MPAS-Ocean runoff fluxes (iceRunoffFlux and
riverRunoffFlux) area-integrated over a Greenland bounding box, comparing
a baseline simulation (no GrIS coupling) against a test simulation with
GrIS SMB coupling enabled.

Data source: *.mpaso.hist.am.timeSeriesStatsMonthly.*.nc

Usage
-----
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    python plot_timeseries_runoff.py

Output is saved to ./plots/runoff_timeseries.png
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

# Ensure utils.py in the same directory is importable.
_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)

import utils

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CTRL_RUN_DIR = (
    '/pscratch/sd/s/sprice/e3sm_scratch/pm-cpu/'
    '20260320.WCYCL2010NS.ne30pg2_r05_IcoswISC30E3r5.pm-cpu.baseline/run/'
)
TEST_RUN_DIR = (
    '/pscratch/sd/s/sprice/e3sm_scratch/pm-cpu/'
    '20260305.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40.pm-cpu.testConfigNewSMBandIC/run/'
)

CTRL_LABEL = 'Baseline'
TEST_LABEL = 'with GIS coupling)'

# Greenland bounding box
LAT_MIN = 55.0   # degrees N
LAT_MAX = 85.0   # degrees N
LON_MIN = -75.0  # degrees (W is negative)
LON_MAX = -15.0  # degrees

# Fallback mesh/init file (if no mpaso.rst found in run dirs)
MESH_FALLBACK = (
    '/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/IcoswISC30E3r5/'
    'mpaso.IcoswISC30E3r5.20231120.nc'
)

OUTPUT_DIR = './plots'

# Year range to process (None = no limit)
YEAR_START = 2011
YEAR_END = 2100

# Conversion factor: kg/s -> Gt/yr
KG_PER_S_TO_GT_PER_YR = 365.25 * 86400.0 / 1e12

# Y-axis limits for each panel: (ymin, ymax).  Use None for auto.
YLIM_RIVER = (0, None)   # riverRunoffFlux panel
YLIM_ICE   = (0, None)   # iceRunoffFlux panel

# Plot styling
COLORS = {'ctrl': 'C0', 'test': 'C1'}
LINEWIDTHS = {'ctrl': 1.0, 'test': 1.0}
LINESTYLES = {'ctrl': '-', 'test': '-'}


# ---------------------------------------------------------------------------
# Mesh loading with fallback
# ---------------------------------------------------------------------------

def load_mesh(run_dir):
    """
    Load MPAS-O mesh coordinates. Try mpaso.rst in run_dir first;
    fall back to the static init file.
    """
    try:
        lat, lon, area = utils.load_mpaso_mesh(run_dir)
        return lat, lon, area
    except RuntimeError:
        pass

    print(f"  No mpaso restart in run dir; using fallback mesh: {MESH_FALLBACK}")
    ds = xr.open_dataset(MESH_FALLBACK, decode_times=False, mask_and_scale=False)

    lat = None
    for v in ['latCell', 'lat']:
        if v in ds:
            lat = utils.rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lat is None:
        raise RuntimeError(f"No latCell in {MESH_FALLBACK}")

    lon = None
    for v in ['lonCell', 'lon']:
        if v in ds:
            lon = utils.rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lon is None:
        raise RuntimeError(f"No lonCell in {MESH_FALLBACK}")
    lon = utils.normalize_lon(lon)

    area = None
    for v in ['areaCell', 'area']:
        if v in ds:
            area = ds[v].values.ravel().astype(float)
            break
    if area is None:
        raise RuntimeError(f"No areaCell in {MESH_FALLBACK}")

    ds.close()
    print(f"    Mesh: {len(lat)} cells (from fallback init file)")
    return lat, lon, area


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_runoff_timeseries(run_dir, lat, lon, area, mask):
    """
    Load MPAS-O timeSeriesStatsMonthly files and compute area-integrated
    iceRunoffFlux and riverRunoffFlux over the masked region.

    Parameters
    ----------
    run_dir : str
        Path to the simulation run directory.
    lat, lon, area : 1-D arrays
        Mesh coordinates and cell areas (from load_mesh).
    mask : 1-D bool array
        True for cells within the Greenland bounding box.

    Returns
    -------
    pd.DataFrame with columns 'iceRunoffFlux_Gt_yr' and 'riverRunoffFlux_Gt_yr',
    indexed by datetime. Returns None if no files found.
    """
    files = utils.find_files(
        run_dir, '*.mpaso.hist.am.timeSeriesStatsMonthly.*.nc'
    )
    if not files:
        return None

    # Filter files by year range (date is in filename: ...YYYY-MM-DD.nc)
    if YEAR_START is not None or YEAR_END is not None:
        import re
        _date_re = re.compile(r'\.(\d{4})-\d{2}-\d{2}\.nc$')
        filtered = []
        for f in files:
            m = _date_re.search(f)
            if m:
                yr = int(m.group(1))
                if YEAR_START is not None and yr < YEAR_START:
                    continue
                if YEAR_END is not None and yr > YEAR_END:
                    continue
            filtered.append(f)
        files = filtered

    if not files:
        print("    No files in requested year range.")
        return None

    area_mask = area[mask]
    n_mask = mask.sum()
    print(f"    Cells in Greenland box: {n_mask} / {len(lat)}")

    all_times = []
    ice_flux = []
    river_flux = []

    print(f"    Processing {len(files)} file(s)...")
    for fpath in files:
        try:
            ds = xr.open_dataset(fpath, decode_times=False, mask_and_scale=True)
        except Exception as exc:
            print(f"    Warning: cannot open {os.path.basename(fpath)}: {exc}")
            continue

        n_t = ds.sizes.get('Time', 0)

        # Parse time
        t = utils.parse_mpas_xtime(ds)
        if t is not None and len(t) == n_t:
            all_times.extend(t)
        else:
            all_times.extend([None] * n_t)

        # Ice runoff flux
        ice_var = None
        for vname in ['timeMonthly_avg_iceRunoffFlux', 'iceRunoffFlux']:
            if vname in ds:
                ice_var = vname
                break
        if ice_var is not None:
            data = ds[ice_var].values[:, mask].astype(float)  # (Time, nMask)
            data = np.where(np.abs(data) < 1e30, data, 0.0)
            integrated = np.sum(data * area_mask[np.newaxis, :], axis=1)
            ice_flux.extend((integrated * KG_PER_S_TO_GT_PER_YR).tolist())
        else:
            ice_flux.extend([0.0] * n_t)

        # River runoff flux
        river_var = None
        for vname in ['timeMonthly_avg_riverRunoffFlux', 'riverRunoffFlux']:
            if vname in ds:
                river_var = vname
                break
        if river_var is not None:
            data = ds[river_var].values[:, mask].astype(float)
            data = np.where(np.abs(data) < 1e30, data, 0.0)
            integrated = np.sum(data * area_mask[np.newaxis, :], axis=1)
            river_flux.extend((integrated * KG_PER_S_TO_GT_PER_YR).tolist())
        else:
            river_flux.extend([0.0] * n_t)

        ds.close()

    # Assemble DataFrame
    df = pd.DataFrame({
        'iceRunoffFlux_Gt_yr': ice_flux,
        'riverRunoffFlux_Gt_yr': river_flux,
    })

    if all(t is not None for t in all_times):
        df.index = pd.DatetimeIndex(all_times)
    else:
        df.index = pd.RangeIndex(len(df))
        print("    Warning: partial time axis; using integer index.")

    print(f"    Loaded {len(df)} monthly records.")
    return df


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_plot(ctrl_df, test_df, output_dir):
    """Create a 2-panel time series plot and save to output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)

    panels = [
        ('riverRunoffFlux_Gt_yr', 'River Runoff Flux (Greenland region)', YLIM_RIVER),
        ('iceRunoffFlux_Gt_yr', 'Ice Runoff Flux (Greenland region)', YLIM_ICE),
    ]

    for ax, (col, title, ylim) in zip(axes, panels):
        if ctrl_df is not None and col in ctrl_df.columns:
            ax.plot(ctrl_df.index, ctrl_df[col],
                    color=COLORS['ctrl'], lw=LINEWIDTHS['ctrl'],
                    ls=LINESTYLES['ctrl'], label=CTRL_LABEL, zorder=2)
        if test_df is not None and col in test_df.columns:
            ax.plot(test_df.index, test_df[col],
                    color=COLORS['test'], lw=LINEWIDTHS['test'],
                    ls=LINESTYLES['test'], label=TEST_LABEL, zorder=3)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel('Flux (Gt/yr)')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        # Apply y-axis limits
        if ylim[0] is not None:
            ax.set_ylim(bottom=ylim[0])
        if ylim[1] is not None:
            ax.set_ylim(top=ylim[1])

    # Format x-axis
    axes[-1].set_xlabel('Year')
    axes[-1].xaxis.set_major_locator(mdates.YearLocator(10))
    axes[-1].xaxis.set_minor_locator(mdates.YearLocator(5))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    fig.suptitle(
        f'MPAS-Ocean Runoff Fluxes — Greenland '
        f'[{LAT_MIN}–{LAT_MAX}°N, {LON_MIN}–{LON_MAX}°E]',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    outpath = os.path.join(output_dir, 'runoff_timeseries.png')
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("MPAS-Ocean Runoff Flux Time Series")
    print("=" * 70)

    # --- Load mesh (use control run dir first, then fallback) ---
    print("\n[1/3] Loading MPAS-O mesh...")
    lat, lon, area = load_mesh(CTRL_RUN_DIR)

    # --- Build Greenland bounding-box mask ---
    mask = (
        (lat >= LAT_MIN) & (lat <= LAT_MAX) &
        (lon >= LON_MIN) & (lon <= LON_MAX)
    )
    print(f"  Greenland mask: {mask.sum()} cells "
          f"(lat [{LAT_MIN}, {LAT_MAX}], lon [{LON_MIN}, {LON_MAX}])")

    if mask.sum() == 0:
        print("  ERROR: No cells in Greenland box. Check coordinate ranges.")
        sys.exit(1)

    # --- Load data from both simulations ---
    print(f"\n[2/3] Loading monthly runoff data...")
    print(f"  Control: {CTRL_RUN_DIR}")
    ctrl_df = load_runoff_timeseries(CTRL_RUN_DIR, lat, lon, area, mask)

    print(f"  Test: {TEST_RUN_DIR}")
    test_df = load_runoff_timeseries(TEST_RUN_DIR, lat, lon, area, mask)

    if ctrl_df is None and test_df is None:
        print("  ERROR: No data loaded from either simulation.")
        sys.exit(1)

    # --- Plot ---
    print(f"\n[3/3] Generating plot...")
    make_plot(ctrl_df, test_df, OUTPUT_DIR)

    print("\nDone.")


if __name__ == '__main__':
    main()
