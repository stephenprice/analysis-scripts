#!/usr/bin/env python3
"""
mali_03_fast_cell_report.py
----------------------------
For each MALI monthly history file spanning the sub-daily-timestep
transition, identifies the cells most likely to be causing the CFL
constraint.  Reports their index, geographic location (lat/lon degrees),
surface speed (m/yr), ice thickness, effective cell spacing, and
estimated advective-CFL allowable dt.

The per-cell CFL estimate mirrors the actual MALI criterion:
    stable_dt(cell) = min over edges e of cell:
                          CFL_FRACTION * dcEdge[e] / |normalVelocity[e]|

    CFL_FRACTION is read live from mali_in
    (config_adaptive_timestep_cfl_fraction).
Cells where stable_dt < 86400 s (1 day) are flagged.

Output per date:
  - Printed rank table (top TOP_N cells, sorted by speed descending)
  - CSV with columns including speed_myr, max_normalvel_myr, stable_dt_days

Final summary printed across all dates.

Run:
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    cd /global/cfs/cdirs/e3sm/sprice/scripts/mali_fast_flow
    python mali_03_fast_cell_report.py
"""

import csv
import glob
import os
import re
import numpy as np
import xarray as xr

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
MALI_IN      = os.path.join(RUN_DIR, "mali_in")

SPY          = 3.1536e7   # s / yr  (m/s → m/yr)
SECS_PER_DAY = 86400.0
ICE_MASK_BIT = 32          # li_mask_ValueIce
TOP_N        = 50          # number of fastest cells to report per date
ALBANY_BIT   = 64          # li_mask_ValueAlbanyActive
FLOATING_BIT = 4           # li_mask_ValueFloating
SIM_START_YR = 2010
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

# ── read model configuration from mali_in ────────────────────────────────────
def _read_namelist_float(filepath, param):
    """Return the first matching float value from a Fortran namelist file."""
    pattern = re.compile(
        r"^\s*" + re.escape(param) + r"\s*=\s*([0-9.eE+\-]+)", re.IGNORECASE
    )
    with open(filepath) as fh:
        for line in fh:
            m = pattern.match(line)
            if m:
                return float(m.group(1))
    raise KeyError(f"'{param}' not found in {filepath}")

CFL_FRACTION = _read_namelist_float(MALI_IN, "config_adaptive_timestep_cfl_fraction")
print(f"CFL fraction (from mali_in): {CFL_FRACTION}")

# ── load mesh topology from a mali restart file ───────────────────────────
rst_files = sorted(glob.glob(os.path.join(RUN_DIR, f"{CASE}.mali.rst.*.nc")))
if not rst_files:
    raise FileNotFoundError(
        f"No mali restart files found in {RUN_DIR}. "
        "Need one for mesh topology (latCell, lonCell, areaCell)."
    )
MESH_FILE = rst_files[0]
print(f"Mesh topology source: {MESH_FILE}")
ds_mesh   = xr.open_dataset(MESH_FILE, decode_cf=False)[
    ["latCell", "lonCell", "areaCell", "edgesOnCell", "nEdgesOnCell", "dcEdge"]
]
lat_deg_all = np.degrees(ds_mesh["latCell"].values.astype(float))
lon_deg_all = np.degrees(ds_mesh["lonCell"].values.astype(float))
lon_deg_all = ((lon_deg_all + 180.0) % 360.0) - 180.0
area_cell   = ds_mesh["areaCell"].values.astype(float)
dc_edge     = ds_mesh["dcEdge"].values.astype(float)
edges_on_cell   = ds_mesh["edgesOnCell"].values.astype(int) - 1  # 0-indexed
n_edges_on_cell = ds_mesh["nEdgesOnCell"].values.astype(int)
_max_edges  = edges_on_cell.shape[1]
edge_mask   = np.arange(_max_edges)[None, :] < n_edges_on_cell[:, None]  # (nCells, maxEdges)
safe_edges  = np.where(edge_mask, edges_on_cell, 0)  # replace padding with edge 0
# Pre-compute per-cell minimum dcEdge for the surfaceSpeed fallback path
_dc_per_cell = dc_edge[safe_edges].copy()           # (nCells, maxEdges)
_dc_per_cell[~edge_mask] = np.inf
min_dc_cell = _dc_per_cell.min(axis=1)              # (nCells,)
ds_mesh.close()
print(f"  nCells = {len(area_cell)}, nEdges = {len(dc_edge)}")

# ── discover available history files and let user choose a window ───────────
hist_files = sorted(glob.glob(
    os.path.join(RUN_DIR, f"{CASE}.mali.hist.????-??-??_00000.nc")
))
all_dates = []
for fp in hist_files:
    bn        = os.path.basename(fp)
    date_part = bn.split(".mali.hist.")[-1].replace("_00000.nc", "")
    yr        = int(date_part[:4])
    mo        = int(date_part[5:7])
    sim_yr    = (yr - SIM_START_YR) + (mo - 1) / 12.0
    all_dates.append((date_part, sim_yr))

print(f"\nFound {len(all_dates)} monthly history files.")
print(f"Simulation year range: {all_dates[0][1]:.2f} – {all_dates[-1][1]:.2f}")
print()

print("Enter a time window (calendar years, YYYY) to analyse.")
print(f"(Simulation started {SIM_START_YR}; available range "
      f"~{int(SIM_START_YR + all_dates[0][1])}–{int(SIM_START_YR + all_dates[-1][1])})")
_lo = input("  Window start calendar year (YYYY): ").strip()
_hi = input("  Window end   calendar year (YYYY): ").strip()
try:
    WIN_LO = float(_lo) - SIM_START_YR
    WIN_HI = float(_hi) - SIM_START_YR + 1.0
except ValueError:
    raise ValueError("Please enter 4-digit calendar years (e.g. 2023).")

DATES = [d for d, y in all_dates if WIN_LO <= y <= WIN_HI]
print(f"\nAnalysing {len(DATES)} files in "
      f"{int(WIN_LO+SIM_START_YR)}–{int(WIN_HI+SIM_START_YR)} "
      f"(sim yr {WIN_LO:.1f}–{WIN_HI:.1f}):")
for d in DATES:
    print(f"  {d}")
print()

summary_rows = []

for date_str in DATES:
    fname = os.path.join(RUN_DIR, f"{CASE}.mali.hist.{date_str}_00000.nc")
    yr_since_start = (int(date_str[:4]) - SIM_START_YR) + int(date_str[5:7]) / 12.0

    if not os.path.exists(fname):
        print(f"File not found, skipping: {date_str}")
        continue

    print(f"\n{'='*72}")
    print(f"  {date_str}  (sim yr {yr_since_start:.2f})")
    print(f"  {fname}")
    print(f"{'='*72}")

    ds = xr.open_dataset(fname, decode_cf=False)

    # -- mesh geometry loaded from restart file (not in history files)
    lat_deg   = lat_deg_all
    lon_deg   = lon_deg_all
    # (dx_eff and area_cell already loaded from mesh)

    # -- cell mask and derived flags (time index 0)
    cell_mask    = ds["cellMask"].values[0, :].astype(int)
    ice_cells    = (cell_mask & ICE_MASK_BIT)   != 0
    albany_cells = (cell_mask & ALBANY_BIT)     != 0
    float_cells  = (cell_mask & FLOATING_BIT)   != 0

    # -- surface speed: m/s (Registry) → m/yr
    speed_mps = ds["surfaceSpeed"].values[0, :].astype(float)
    speed_myr = speed_mps * SPY

    # -- normalVelocity and per-cell stable CFL dt.
    # Use normalVelocity at edges if available; fall back to surfaceSpeed.
    if "normalVelocity" in ds:
        norm_vel_raw = ds["normalVelocity"].values[0]  # (nEdges,) or (nEdges, nVertLevels)
        if norm_vel_raw.ndim == 2:
            nv_abs = np.abs(norm_vel_raw).max(axis=-1)
        else:
            nv_abs = np.abs(norm_vel_raw)
        nv_per_cell_abs = nv_abs[safe_edges]           # (nCells, maxEdges)
        nv_per_cell_abs[~edge_mask] = 0.0
        max_nv_mps = nv_per_cell_abs.max(axis=1)       # (nCells,)
        dc_per = dc_edge[safe_edges]
        with np.errstate(divide="ignore", invalid="ignore"):
            dt_per = np.where(
                edge_mask & (nv_per_cell_abs > 0),
                CFL_FRACTION * dc_per / nv_per_cell_abs,
                1e15,
            )
        dt_acfl_est = dt_per.min(axis=1)
    else:
        # normalVelocity not in history file: use surfaceSpeed with min dcEdge
        # (conservative: speed magnitude >= edge normal component)
        max_nv_mps = speed_mps.copy()                  # best available proxy
        with np.errstate(divide="ignore", invalid="ignore"):
            dt_acfl_est = np.where(
                speed_mps > 0,
                CFL_FRACTION * min_dc_cell / speed_mps,
                1e15,
            )
    max_nv_myr = max_nv_mps * SPY

    # -- thickness
    thickness = ds["thickness"].values[0, :].astype(float)

    # -- global scalars at this output time
    deltat_val  = float(ds["deltat"].values.ravel()[0])
    acfl_global = float(ds["allowableDtACFL"].values.ravel()[0])
    dcfl_global = float(ds["allowableDtDCFL"].values.ravel()[0])

    print(f"  Global deltat        : {deltat_val:.0f} s  "
          f"({deltat_val/3600:.4f} h)")
    print(f"  allowableDtACFL      : {acfl_global:.0f} s  "
          f"({acfl_global/3600:.4f} h)")
    print(f"  allowableDtDCFL      : {dcfl_global:.0f} s  "
          f"({dcfl_global/3600:.4f} h)")

    # -- restrict to ice-bearing cells
    ice_idx = np.where(ice_cells)[0]
    if len(ice_idx) == 0:
        print("  No ice cells found.")
        ds.close()
        continue

    speed_ice   = speed_myr[ice_idx]
    dt_ice      = dt_acfl_est[ice_idx]
    nv_ice      = max_nv_myr[ice_idx]
    lat_ice     = lat_deg[ice_idx]
    lon_ice     = lon_deg[ice_idx]
    thick_ice   = thickness[ice_idx]
    dx_ice      = dc_edge[safe_edges][ice_idx].max(axis=1)  # representative edge length
    alb_ice     = albany_cells[ice_idx]
    flt_ice     = float_cells[ice_idx]
    mask_ice    = cell_mask[ice_idx]

    n_ice       = len(ice_idx)
    n_cfl       = int((dt_ice < SECS_PER_DAY).sum())
    n_albany    = int(alb_ice.sum())

    print(f"  Ice cells            : {n_ice}")
    print(f"  Albany-active cells  : {n_albany}")
    print(f"  CFL-limited (est.)   : {n_cfl}  (stable dt < 1 day, frac={CFL_FRACTION})")
    print(f"  Max speed            : {speed_ice.max():.4e} m/yr")
    print(f"  Max |normalVelocity| : {nv_ice.max():.4e} m/yr")
    print(f"  Min est. CFL dt      : {dt_ice.min()/3600:.4f} h")
    print()

    # Sort by speed descending → smallest CFL dt first
    order   = np.argsort(speed_ice)[::-1]
    n_print = min(TOP_N, len(order))

    # Print table
    col_w = 65
    hdr = (f"  {'rnk':>4}  {'cellIdx':>8}  {'lat°':>7}  {'lon°':>8}  "
           f"{'spd_myr':>12}  {'maxNV_myr':>12}  {'dt_days':>9}  {'thk_m':>7}  "
           f"{'dx_km':>6}  {'Albany':>6}  {'Float':>5}  CFL?")
    sep  = "  " + "-" * (len(hdr) - 2)
    print(hdr)
    print(sep)

    csv_rows = []
    for rank in range(n_print):
        j         = order[rank]
        cell_glob = ice_idx[j]
        spd       = speed_ice[j]
        nv        = nv_ice[j]
        dt_days   = dt_ice[j] / SECS_PER_DAY
        lat       = lat_ice[j]
        lon       = lon_ice[j]
        thk       = thick_ice[j]
        dx        = dx_ice[j]
        alb_flag  = "YES" if alb_ice[j] else "no"
        flt_flag  = "YES" if flt_ice[j] else "no"
        cfl_flag  = "YES" if dt_ice[j] < SECS_PER_DAY else "no"

        print(
            f"  {rank+1:>4}  {cell_glob:>8d}  {lat:>7.3f}  {lon:>8.3f}  "
            f"{spd:>12.4e}  {nv:>12.4e}  {dt_days:>9.4f}  {thk:>7.1f}  "
            f"{dx/1000:>6.2f}  {alb_flag:>6}  {flt_flag:>5}  {cfl_flag}"
        )

        csv_rows.append({
            "rank":              rank + 1,
            "cellIndex":         int(cell_glob),
            "lat_deg":           round(float(lat), 5),
            "lon_deg":           round(float(lon), 5),
            "speed_myr":         float(spd),
            "max_normalvel_myr": float(nv),
            "stable_dt_h":       round(float(dt_ice[j] / 3600.0), 6),
            "stable_dt_days":    round(float(dt_days), 6),
            "thickness_m":       round(float(thk), 2),
            "dx_eff_km":         round(float(dx / 1000), 3),
            "albanyActive":      alb_flag,
            "floating":          flt_flag,
            "cfl_limited":       cfl_flag,
            "cellMask_raw":      int(mask_ice[j]),
        })

    # Save CSV
    csv_file = os.path.join(OUT_DIR, f"mali_fast_cells_{date_str}.csv")
    with open(csv_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_rows[0].keys())
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"\n  CSV written: {csv_file}")

    summary_rows.append((
        date_str,
        f"{yr_since_start:.2f}",
        f"{deltat_val/3600:.4f}",
        f"{acfl_global/3600:.4f}",
        n_cfl,
        f"{speed_ice.max():.4e}",
        f"{dt_ice.min()/3600:.4f}",
    ))
    ds.close()

# ── final summary ─────────────────────────────────────────────────────────────
print(f"\n\n{'='*72}")
print("  SUMMARY — all dates")
print(f"{'='*72}")
hdr = (f"  {'date':12}  {'sim_yr':>7}  {'dt_h':>8}  "
       f"{'aCFL_h':>8}  {'n_CFL':>7}  {'max_spd_myr':>13}  {'min_dtEst_h':>12}")
print(hdr)
print("  " + "-" * (len(hdr) - 2))
for r in summary_rows:
    print(f"  {r[0]:12}  {r[1]:>7}  {r[2]:>8}  "
          f"{r[3]:>8}  {r[4]:>7}  {r[5]:>13}  {r[6]:>12}")
print()
print(f"Note: n_CFL = ice cells where stable dt (= dx*{CFL_FRACTION}/speed) < 1 day.")
print("      This uses the same CFL fraction as config_adaptive_timestep_cfl_fraction.")
print("      The actual limiting cells are those that set allowableDtACFL in the solver.")
