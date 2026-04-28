#!/usr/bin/env python3
"""
mali_01_globalStats_timeseries.py
----------------------------------
Plot MALI deltat and surfaceSpeedMax time series from the analysis-member
globalStats file.  Identifies the exact first occurrence of a sub-daily
timestep (deltat < 86400 s) and prints a detailed table centred on that
event.

Run:
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    cd /global/cfs/cdirs/e3sm/sprice/scripts/mali_fast_flow
    python mali_01_globalStats_timeseries.py
"""

import os
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── configuration ─────────────────────────────────────────────────────────────
RUN_DIR = (
    "/pscratch/sd/s/sprice/e3sm_scratch/pm-cpu/"
    "20260305.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40.pm-cpu.testConfigNewSMBandIC/run"
)
CASE = (
    "20260305.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40."
    "pm-cpu.testConfigNewSMBandIC"
)
OUT_DIR      = "/global/cfs/cdirs/e3sm/sprice/analysis/mali_fast_flow"
SIM_START_YR = 2010
SECS_PER_DAY = 86400.0
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

# ── user-specified time window for transition search ──────────────────────────
print("Time window for sub-daily transition search")
print(f"  Enter calendar years YYYY (simulation started {SIM_START_YR}).")
print("  Leave blank to search the entire record.")
_lo = input("  Start calendar year (YYYY) [default: start]: ").strip()
_hi = input("  End   calendar year (YYYY) [default: end  ]: ").strip()
SEARCH_YR_LO = (float(_lo) - SIM_START_YR) if _lo else 0.0
SEARCH_YR_HI = (float(_hi) - SIM_START_YR + 1.0) if _hi else np.inf
_hi_str = str(int(SIM_START_YR + SEARCH_YR_HI)) if np.isfinite(SEARCH_YR_HI) else 'end'
print(f"  Searching {int(SIM_START_YR + SEARCH_YR_LO)} – {_hi_str}"
      f" (sim yr {SEARCH_YR_LO:.1f} – "
      + (f"{SEARCH_YR_HI:.1f})" if np.isfinite(SEARCH_YR_HI) else "end)"))
print()

stats_file = os.path.join(
    RUN_DIR, f"{CASE}.mali.hist.am.globalStats.0001-01-01_00000.nc"
)
print(f"Reading: {stats_file}")
ds = xr.open_dataset(stats_file, decode_cf=False)

print("Variables in file:")
for vn, v in ds.data_vars.items():
    print(f"  {vn:40s}  dims={v.dims}  shape={v.shape}")
print()

# ── time axis ─────────────────────────────────────────────────────────────────
days      = ds["daysSinceStart"].values.astype(float).ravel()
years_sim = days / 365.0

# ── deltat (seconds per mali step) ───────────────────────────────────────────
deltat = ds["deltat"].values.astype(float).ravel()

# ── surfaceSpeedMax ───────────────────────────────────────────────────────────
# This is written by the globalStats analysis member in m/yr
# (domain-wide maximum; capped at 5e5 m/yr)
spd_var     = ds["surfaceSpeedMax"]
spd_units   = getattr(spd_var, "units", "m yr^-1")
speedmax    = spd_var.values.astype(float).ravel()
print(f"surfaceSpeedMax units attribute: '{spd_units}'")

# ── locate first sub-daily timestep within the user window ──────────────────
window_mask = (years_sim >= SEARCH_YR_LO) & (years_sim <= SEARCH_YR_HI)
sub_idx = np.where((deltat < SECS_PER_DAY) & window_mask)[0]

print(f"\n{'='*64}")
print("  TRANSITION: first sub-daily MALI timestep (deltat < 86400 s)")
print(f"{'='*64}")
if len(sub_idx) > 0:
    i0     = sub_idx[0]
    yr0    = years_sim[i0]
    cal_yr = SIM_START_YR + int(yr0)
    cal_mo = int((yr0 % 1) * 12) + 1
    print(f"  Index           : {i0}")
    print(f"  Simulation day  : {days[i0]:.1f}")
    print(f"  Simulation year : {yr0:.4f}")
    print(f"  Approx date     : {cal_yr}-{cal_mo:02d}")
    print(f"  deltat          : {deltat[i0]:.0f} s  "
          f"({deltat[i0]/3600:.4f} h)")
    print(f"  surfaceSpeedMax : {speedmax[i0]:.6e} {spd_units}")
    print()

    # Table: 10 rows before and 20 rows after transition
    lo = max(0, i0 - 10)
    hi = min(len(days), i0 + 21)
    hdr = (f"  {'idx':>7}  {'sim_day':>10}  {'sim_yr':>8}  "
           f"{'deltat_h':>10}  {'speedMax':>14}  flag")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for i in range(lo, hi):
        flag = "  ◄ FIRST SUB-DAILY" if i == i0 else ""
        print(f"  {i:>7d}  {days[i]:>10.1f}  {years_sim[i]:>8.4f}  "
              f"{deltat[i]/3600:>10.5f}  {speedmax[i]:>14.6e}{flag}")
else:
    print("  No sub-daily timesteps found in this dataset.")
print()

# ── plot ───────────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

# -- deltat panel
ax1.plot(years_sim, deltat / 3600, color="steelblue", lw=0.6, alpha=0.8)
ax1.axhline(24, color="crimson", ls="--", lw=1.5, label="24 h (1 day)")
ax1.set_yscale("log")
ax1.set_ylabel("deltat  [hours]", fontsize=11)
ax1.set_title("MALI adaptive timestep", fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=9, loc="upper right")
if len(sub_idx):
    ax1.axvline(years_sim[sub_idx[0]], color="darkorange", ls=":", lw=2,
                label=f"First sub-daily  (sim yr {years_sim[sub_idx[0]]:.3f})")
    ax1.legend(fontsize=9, loc="upper right")

# -- surfaceSpeedMax panel
ax2.plot(years_sim, speedmax, color="firebrick", lw=0.6, alpha=0.8)
ax2.set_yscale("log")
ax2.set_ylabel(f"surfaceSpeedMax  [{spd_units}]", fontsize=11)
ax2.set_xlabel("Simulation years elapsed", fontsize=11)
ax2.set_title("MALI domain-wide maximum surface speed", fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=9, loc="upper right")
if len(sub_idx):
    ax2.axvline(years_sim[sub_idx[0]], color="darkorange", ls=":", lw=2,
                label=f"First sub-daily  (sim yr {years_sim[sub_idx[0]]:.3f})")
    ax2.legend(fontsize=9, loc="upper right")

# Dual-label x-axis: sim year + approx calendar year
# Zoom to user window (with a little margin)
xlo = max(0, SEARCH_YR_LO - 1)
xhi = min(years_sim[-1], SEARCH_YR_HI + 1) if np.isfinite(SEARCH_YR_HI) else years_sim[-1]
ax1.set_xlim(xlo, xhi)
ax2.set_xlim(xlo, xhi)

tick_step = max(1, int((xhi - xlo) / 10))
tick_yrs = np.arange(int(xlo), int(xhi) + tick_step, tick_step).astype(float)
ax2.set_xticks(tick_yrs)
ax2.set_xticklabels(
    [f"{y:.0f}\n({int(SIM_START_YR + y)})" for y in tick_yrs], fontsize=8
)

plt.suptitle(
    f"MALI Global Statistics — {CASE}",
    fontsize=8, y=1.01
)
plt.tight_layout()

outfile = os.path.join(OUT_DIR, "mali_01_globalStats_timeseries.png")
plt.savefig(outfile, dpi=150, bbox_inches="tight")
print(f"Plot saved: {outfile}")
ds.close()
