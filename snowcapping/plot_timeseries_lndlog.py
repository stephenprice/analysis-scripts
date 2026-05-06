"""
plot_timeseries_lndlog.py
=========================
Plot SNOWCAP_LATBAND_DIAG and SNOWCAP_LATBAND_DIAG_CUM diagnostics from the
ELM land-model log files, comparing the fix run against the control run.

Two figures are saved:
  lndlog_step_diagnostics.png — per-diagnostic-step values
  lndlog_cum_diagnostics.png  — cumulative values

Each figure has five panels:
  ice_mass_kg       — ice mass removed by snowcapping (kg)
  latent_energy_j   — latent energy removed (J)
  mean_cooling_k    — area-weighted mean column cooling (K)
  max_cooling_k     — maximum single-column cooling (K)
  cols              — number of affected land columns

The fix run shows its actual logged values.
The control run (where the diagnostic is inactive) is plotted as a flat
zero baseline on the same axes so the difference is immediately apparent.

Data source: lnd.log.*.gz files in each run directory.

Usage
-----
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_login.sh
    python plot_timeseries_lndlog.py

Edit config.py to change run directories, labels, or styling.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)

import config
import utils


# ---------------------------------------------------------------------------
# Column metadata for the two figures
# ---------------------------------------------------------------------------

_STEP_COLS = [
    'step_ice_mass_kg',
    'step_latent_energy_j',
    'step_mean_cooling_k',
    'step_max_cooling_k',
    'step_cols',
]
_CUM_COLS = [
    'cum_ice_mass_kg',
    'cum_latent_energy_j',
    'cum_mean_cooling_k',
    'cum_max_cooling_k',
    'cum_cols',
]

_STEP_LABELS = {
    'step_ice_mass_kg':      'Ice Mass Removed (kg)',
    'step_latent_energy_j':  'Latent Energy Removed (J)',
    'step_mean_cooling_k':   'Area-wtd Mean Cooling (K)',
    'step_max_cooling_k':    'Max Column Cooling (K)',
    'step_cols':             'Affected Land Columns (count)',
}
_CUM_LABELS = {
    'cum_ice_mass_kg':      'Cumulative Ice Mass Removed (kg)',
    'cum_latent_energy_j':  'Cumulative Latent Energy Removed (J)',
    'cum_mean_cooling_k':   'Area-wtd Mean Cooling — cumulative avg (K)',
    'cum_max_cooling_k':    'Running Max Column Cooling (K)',
    'cum_cols':             'Cumulative Affected Land Columns (count)',
}

_TITLES = {
    'step': 'SNOWCAP_LATBAND_DIAG  (per diagnostic step)',
    'cum':  'SNOWCAP_LATBAND_DIAG_CUM  (cumulative)',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _zero_df_like(ref_df, cols):
    """Return a DataFrame of zeros indexed identically to ref_df."""
    return pd.DataFrame(0.0, index=ref_df.index, columns=cols)


def _has_datetime_index(df):
    return df is not None and not df.empty and isinstance(df.index, pd.DatetimeIndex)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_lndlog(step_fix, step_ctrl, cum_fix, cum_ctrl, output_dir):
    """
    Generate and save the two lndlog diagnostic figures.

    Parameters
    ----------
    step_fix, step_ctrl  : DataFrames from parse_lnd_logs (per-step entries)
    cum_fix,  cum_ctrl   : DataFrames from parse_lnd_logs (cumulative entries)
    output_dir           : directory in which to save the PNGs
    """
    os.makedirs(output_dir, exist_ok=True)

    specs = [
        ('step', 'lndlog_step_diagnostics', _STEP_COLS, _STEP_LABELS,
         step_fix, step_ctrl),
        ('cum',  'lndlog_cum_diagnostics',  _CUM_COLS,  _CUM_LABELS,
         cum_fix,  cum_ctrl),
    ]

    for kind, fig_stem, cols, labels, fix_df, ctrl_df in specs:
        # Determine which cols actually appear in the data
        avail = set()
        if fix_df  is not None and not fix_df.empty:  avail |= set(fix_df.columns)
        if ctrl_df is not None and not ctrl_df.empty: avail |= set(ctrl_df.columns)
        plot_cols = [c for c in cols if c in avail]

        if not plot_cols:
            print(f"  No data for {fig_stem}; skipping.")
            continue

        nrows = len(plot_cols)
        fig, axes = plt.subplots(nrows, 1,
                                 figsize=(13, 2.8 * nrows),
                                 sharex=True)
        if nrows == 1:
            axes = [axes]

        for ax, col in zip(axes, plot_cols):
            label_str = labels.get(col, col)

            # Fix run
            if fix_df is not None and not fix_df.empty and col in fix_df.columns:
                ax.plot(fix_df.index, fix_df[col],
                        color=config.COLORS['fix'],
                        ls=config.LINESTYLES['fix'],
                        lw=config.LINEWIDTHS['fix'],
                        label=config.FIX_LABEL, zorder=3)

            # Control run
            if ctrl_df is not None and not ctrl_df.empty and col in ctrl_df.columns:
                ax.plot(ctrl_df.index, ctrl_df[col],
                        color=config.COLORS['ctrl'],
                        ls=config.LINESTYLES['ctrl'],
                        lw=config.LINEWIDTHS['ctrl'],
                        label=config.CTRL_LABEL, zorder=2)

            ax.set_ylabel(label_str, fontsize=8)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3, linewidth=0.5)

        # Format x-axis
        use_dates = _has_datetime_index(fix_df) or _has_datetime_index(ctrl_df)
        if use_dates:
            axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
            axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
            plt.setp(axes[-1].get_xticklabels(), rotation=30, ha='right')

        axes[-1].set_xlabel('Model Date')
        fig.suptitle(_TITLES[kind], fontsize=12, y=1.01)
        fig.tight_layout()

        outpath = os.path.join(output_dir, f'{fig_stem}.png')
        fig.savefig(outpath, dpi=config.DPI, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("lnd.log Diagnostic Time Series")
    print("=" * 60)
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    # Read simulation metadata from the fix run's namelists
    dtime      = utils.get_dtime(config.FIX_RUN_DIR)
    start_date = utils.get_sim_start_date(config.FIX_RUN_DIR)
    print(f"  Simulation start: {start_date},  dtime: {dtime} s  "
          f"(read from fix run lnd_in / drv_in)")

    print(f"\nParsing fix run logs:  {config.FIX_RUN_DIR}")
    step_fix, cum_fix = utils.parse_lnd_logs(
        config.FIX_RUN_DIR, start_date, dtime)

    print(f"\nParsing ctrl run logs: {config.CTRL_RUN_DIR}")
    step_ctrl, cum_ctrl = utils.parse_lnd_logs(
        config.CTRL_RUN_DIR, start_date, dtime)

    # If the control run produced no diagnostics (expected), substitute a
    # zero-valued DataFrame aligned to the fix run's time index so both
    # simulations appear on the same axes with the control as a flat zero line.
    if step_ctrl.empty and not step_fix.empty:
        step_ctrl = _zero_df_like(step_fix, _STEP_COLS)
        print("  Control run: no SNOWCAP_LATBAND_DIAG entries — "
              "plotting control as zeros.")
    if cum_ctrl.empty and not cum_fix.empty:
        cum_ctrl = _zero_df_like(cum_fix, _CUM_COLS)

    print("\nGenerating plots…")
    plot_lndlog(step_fix, step_ctrl, cum_fix, cum_ctrl, config.OUTPUT_DIR)
    print("Done.\n")


if __name__ == '__main__':
    main()
