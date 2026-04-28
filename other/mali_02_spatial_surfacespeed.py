#!/usr/bin/env python3
"""
mali_02_spatial_surfacespeed.py
--------------------------------
Six-panel mosaic map of MALI surfaceSpeed on the native GIS mesh for
dates bracketing the sub-daily-timestep transition.

Mesh topology (latCell, lonCell, areaCell, verticesOnCell, etc.) is
loaded from a mali restart file, because these variables are NOT written
to the monthly history files.

surfaceSpeed in the monthly history files is in m/s (MALI Registry).
It is converted to m/yr for display (multiply by SPY = 3.1536e7 s/yr).

Cells whose stable advective-CFL timestep is < 86400 s (1 day)
are overlaid with a gold fill and black border so their location is clearly visible.

    stable_dt(cell) = min over edges e of cell:
                          CFL_FRACTION * dcEdge[e] / |normalVelocity[e]|

This mirrors the actual MALI criterion: the model uses normalVelocity at cell
edges (not surfaceSpeed) and the CFL fraction from mali_in.

latCell / lonCell in MPAS files are in radians; mosaic.Descriptor
converts them to degrees internally when use_latlon=True.

Run:
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    cd /global/cfs/cdirs/e3sm/sprice/scripts/mali_fast_flow
    python mali_02_spatial_surfacespeed.py
"""

import glob
import os
import re
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import mosaic

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
SPY          = 3.1536e7    # seconds per year (m/s → m/yr)
SECS_PER_DAY = 86400.0
ICE_MASK_BIT = 32          # li_mask_ValueIce from mpas_li_mask.F
SPEED_VMIN   = 10.0        # m/yr  (log scale lower bound)
SPEED_VMAX   = 5e4         # m/yr  (log scale upper bound)
SIM_START_YR = 2010
CFL_BORDER   = (0.0, 1.0, 0.0)  # bright green border for CFL-limited cells (RGB 0,255,0)
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

# ── locate a mali restart file for mesh topology ──────────────────────────────
rst_files = sorted(glob.glob(os.path.join(RUN_DIR, f"{CASE}.mali.rst.*.nc")))
if not rst_files:
    raise FileNotFoundError(
        f"No mali restart files found in {RUN_DIR}. "
        "Need one for mesh topology (latCell, areaCell, etc.)."
    )
MESH_FILE = rst_files[0]
print(f"Mesh topology source: {MESH_FILE}")

# Open the full restart dataset for mosaic — it needs several coordinate and
# connectivity arrays internally and it is simplest to let it find them itself.
# We also extract areaCell here for the per-cell CFL estimate.
ds_mesh = xr.open_dataset(MESH_FILE, decode_cf=False)
# mosaic requires an 'is_periodic' global attribute; MALI restart files
# use 'on_a_sphere' instead.  This is a regional (non-periodic) mesh.
ds_mesh.attrs["is_periodic"] = "NO"
area_cell = ds_mesh["areaCell"].values.astype(float)
# Edge connectivity for per-cell CFL using normalVelocity
edges_on_cell   = ds_mesh["edgesOnCell"].values.astype(int) - 1  # 0-indexed (-1 = padding)
n_edges_on_cell = ds_mesh["nEdgesOnCell"].values.astype(int)
dc_edge         = ds_mesh["dcEdge"].values.astype(float)
_max_edges  = edges_on_cell.shape[1]
edge_mask   = np.arange(_max_edges)[None, :] < n_edges_on_cell[:, None]  # (nCells, maxEdges)
safe_edges  = np.where(edge_mask, edges_on_cell, 0)  # replace padding with edge 0
# Pre-compute per-cell minimum dcEdge for the surfaceSpeed fallback path
dc_per_cell = dc_edge[safe_edges]                    # (nCells, maxEdges)
dc_per_cell[~edge_mask] = np.inf
min_dc_cell = dc_per_cell.min(axis=1)                # (nCells,)
print(f"  Mesh loaded: nCells={len(area_cell)}, nEdges={len(dc_edge)}")

# ── discover available history files and let user choose a window ─────────────
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
    all_dates.append((date_part, sim_yr, fp))

print(f"\nFound {len(all_dates)} monthly history files.")
print(f"Simulation year range: {all_dates[0][1]:.2f} – {all_dates[-1][1]:.2f}")
print()

print("Enter a time window (calendar years, YYYY) to plot.")
print("Six panels will be selected evenly across this window.")
print(f"(Simulation started {SIM_START_YR}; available range "
      f"~{int(SIM_START_YR + all_dates[0][1])}–{int(SIM_START_YR + all_dates[-1][1])})")
_lo = input("  Window start calendar year (YYYY): ").strip()
_hi = input("  Window end   calendar year (YYYY): ").strip()
try:
    WIN_LO = float(_lo) - SIM_START_YR
    WIN_HI = float(_hi) - SIM_START_YR + 1.0
except ValueError:
    raise ValueError("Please enter 4-digit calendar years (e.g. 2023).")

window_dates = [(d, y, fp) for d, y, fp in all_dates if WIN_LO <= y <= WIN_HI]
if len(window_dates) < 2:
    raise ValueError(
        f"Only {len(window_dates)} files found in calendar yr "
        f"[{int(WIN_LO+SIM_START_YR)}, {int(WIN_HI+SIM_START_YR)}]. "
        "Broaden the window."
    )

# Choose up to 6 evenly-spaced dates across the window
n_window = len(window_dates)
if n_window <= 6:
    chosen = window_dates
else:
    idxs   = np.round(np.linspace(0, n_window - 1, 6)).astype(int)
    chosen = [window_dates[i] for i in idxs]

DATES   = [d for d, y, fp in chosen]
SIM_YRS = [y for d, y, fp in chosen]
print(f"\nSelected {len(DATES)} panels:")
for d, y in zip(DATES, SIM_YRS):
    print(f"  {d}  (sim yr {y:.2f})")
print()

# ── map setup ─────────────────────────────────────────────────────────────────
projection = ccrs.NorthPolarStereo(central_longitude=-42.0)
transform  = ccrs.Geodetic()   # mosaic converts radians→degrees internally

norm = mcolors.LogNorm(vmin=SPEED_VMIN, vmax=SPEED_VMAX)
cmap = plt.cm.plasma.copy()
cmap.set_bad(color="none")   # NaN (non-ice) cells → transparent

# ── build mosaic Descriptor from restart file mesh (done once) ─────────────
print("Building mosaic Descriptor from mesh file...")
descriptor = mosaic.Descriptor(ds_mesh, projection, transform, use_latlon=True)
print("  Descriptor ready.")
print()

# ── figure ────────────────────────────────────────────────────────────────────
n_panels = len(DATES)
ncols    = 3
nrows    = int(np.ceil(n_panels / ncols))
fig, axes = plt.subplots(
    nrows, ncols,
    figsize=(6 * ncols, 6 * nrows),
    subplot_kw={"projection": projection},
    constrained_layout=True,
)
axes_flat = np.array(axes).flatten()

for ax in axes_flat[n_panels:]:
    ax.set_visible(False)

for ax, date_str, sim_yr in zip(axes_flat, DATES, SIM_YRS):
    fname = os.path.join(RUN_DIR, f"{CASE}.mali.hist.{date_str}_00000.nc")

    if not os.path.exists(fname):
        ax.set_title(f"{date_str}\n(file not found)", fontsize=9)
        print(f"  Missing: {fname}")
        continue

    print(f"Loading {fname}")
    ds = xr.open_dataset(fname, decode_cf=False)

    # Ice mask
    cell_mask = ds["cellMask"].values[0, :].astype(int)
    ice_cells  = (cell_mask & ICE_MASK_BIT) != 0

    # surfaceSpeed m/s → m/yr
    speed_mps = ds["surfaceSpeed"].values[0, :].astype(float)
    speed_myr = speed_mps * SPY

    # Per-cell stable advective-CFL dt (s).
    # Use normalVelocity at edges if available; fall back to surfaceSpeed.
    if "normalVelocity" in ds:
        norm_vel_raw = ds["normalVelocity"].values[0]  # (nEdges,) or (nEdges, nVertLevels)
        if norm_vel_raw.ndim == 2:
            nv_abs = np.abs(norm_vel_raw).max(axis=-1)
        else:
            nv_abs = np.abs(norm_vel_raw)
        nv_per = nv_abs[safe_edges]                    # (nCells, maxEdges)
        dc_per = dc_edge[safe_edges]
        with np.errstate(divide="ignore", invalid="ignore"):
            dt_per = np.where(
                edge_mask & (nv_per > 0),
                CFL_FRACTION * dc_per / nv_per,
                np.inf,
            )
        dt_acfl_est = dt_per.min(axis=1)
    else:
        # normalVelocity not in history file: use surfaceSpeed with min dcEdge
        # (conservative: speed magnitude >= edge normal component)
        with np.errstate(divide="ignore", invalid="ignore"):
            dt_acfl_est = np.where(
                speed_mps > 0,
                CFL_FRACTION * min_dc_cell / speed_mps,
                np.inf,
            )

    # Global scalars
    acfl_global = float(ds["allowableDtACFL"].values.ravel()[0])
    dcfl_global = float(ds["allowableDtDCFL"].values.ravel()[0])
    deltat_val  = float(ds["deltat"].values.ravel()[0])

    # Mask to ice only for plotting.
    # mosaic.polypcolor requires xarray DataArrays (not bare numpy arrays).
    speed_plot  = xr.DataArray(
        np.where(ice_cells, speed_myr, np.nan), dims=["nCells"]
    )
    cfl_limited = ice_cells & (dt_acfl_est < SECS_PER_DAY)
    n_cfl       = int(cfl_limited.sum())

    # Draw speed field
    mosaic.polypcolor(
        ax, descriptor, speed_plot,
        norm=norm, cmap=cmap, antialiaseds=False,
    )

    # Overlay CFL-limited cells: native speed fill is shown as-is, bright green
    # border only on violating cells. We pass a dummy uniform array so mosaic
    # draws polygons for all cells, set facecolor to none to keep the underlying
    # speed colors visible, then set per-cell edge RGBA so only CFL cells get
    # the green border.
    if n_cfl > 0:
        dummy = xr.DataArray(np.ones(len(cfl_limited)), dims=["nCells"])
        pc_cfl = mosaic.polypcolor(
            ax, descriptor, dummy,
            norm=mcolors.Normalize(vmin=0, vmax=2),
            cmap=mcolors.ListedColormap(["none"]),
            edgecolors="none", linewidths=2.0,
            antialiaseds=True, zorder=6,
        )
        # Per-cell edge RGBA: bright green on CFL cells, fully transparent elsewhere
        edge_rgba = np.zeros((len(cfl_limited), 4), dtype=float)
        edge_rgba[cfl_limited] = (*CFL_BORDER, 1.0)
        pc_cfl.set_edgecolors(edge_rgba)

    ax.add_feature(cfeature.LAND,      facecolor="#d8d8d8", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor="black")
    ax.gridlines(draw_labels=False, linewidth=0.3, color="grey", alpha=0.5)
    ax.set_extent([-70, -10, 58, 86], crs=ccrs.PlateCarree())

    sub_flag = f"  dt={deltat_val/3600:.2f} h" if deltat_val < SECS_PER_DAY else ""
    ax.set_title(
        f"{date_str}  (sim yr {sim_yr:.2f})\n"
        f"CFL-est cells: {n_cfl}{sub_flag}",
        fontsize=9,
    )
    print(f"  dt={deltat_val:.0f}s  aCFL={acfl_global:.0f}s  "
          f"dCFL={dcfl_global:.0f}s  CFL-est cells={n_cfl}")
    ds.close()

# ── shared colorbar and legend ─────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(
    sm, ax=list(axes_flat[:n_panels]), orientation="vertical",
    fraction=0.018, pad=0.02, shrink=0.75,
)
cbar.set_label("surfaceSpeed  [m yr⁻¹]", fontsize=10)

cfl_patch = mpatches.Patch(
    facecolor="none", edgecolor=CFL_BORDER, linewidth=2.0,
    label=f"CFL-limited: stable dt < 1 day  (frac={CFL_FRACTION})",
)
fig.legend(handles=[cfl_patch], loc="lower center", fontsize=9,
           bbox_to_anchor=(0.5, -0.01))

fig.suptitle(
    f"MALI surface speed — {int(WIN_LO+SIM_START_YR)}–{int(WIN_HI+SIM_START_YR)}"
    f" (sim yr {WIN_LO:.1f}–{WIN_HI:.1f})\n{CASE}",
    fontsize=9,
)

outfile = os.path.join(OUT_DIR,
    f"mali_02_surfacespeed_{int(WIN_LO+SIM_START_YR)}-{int(WIN_HI+SIM_START_YR)}.png")
plt.savefig(outfile, dpi=150, bbox_inches="tight")
print(f"\nPlot saved: {outfile}")
ds_mesh.close()
