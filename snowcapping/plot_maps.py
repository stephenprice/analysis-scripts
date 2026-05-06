"""
plot_maps.py
============
Generate monthly map-view plots for five variables, comparing a snowcapping-fix
E3SM run against a control run.

Variables plotted
-----------------
  iceAreaCell         — monthly sea-ice area fraction
                        (mpassi timeSeriesStatsMonthly)
  iceRunoffFlux       — monthly ice runoff flux (kg m⁻² s⁻¹)
                        (mpaso  timeSeriesStatsMonthly)
  riverRunoffFlux     — monthly river runoff flux (kg m⁻² s⁻¹)
                        (mpaso  timeSeriesStatsMonthly)
  QSNWCPICE           — snow-cap ice flux (mm s⁻¹)
                        (elm h0 monthly history)
  SNO_T_davg          — snow temperature depth-averaged over all snow layers (K)
                        (elm h2 monthly history)

Output layout
-------------
One PNG per variable per month: 2 rows × 1 column (control on top, fix on bottom).
For SNO_T_davg, a third row shows the difference (fix − control).
Colourscale is shared across all months and both simulations for each variable
(2nd–98th percentile of all non-NaN values, except iceAreaCell which is [0, 1]).
Colourbars are placed in dedicated space below the map panels (never overlapping).

Files are saved to:
    <OUTPUT_DIR>/maps/<var_name>/month_<MM>_<var_name>.png

  e.g.  plots/maps/iceAreaCell/month_01_iceAreaCell.png
                                month_02_iceAreaCell.png
                                …

Rendering approach
------------------
- MPAS unstructured-mesh variables (iceAreaCell, iceRunoffFlux, riverRunoffFlux)
  are rendered using matplotlib's tripcolor with a masked Delaunay triangulation.
  Triangles whose edges exceed 2× the median cell spacing are masked out, which
  naturally prevents data from being drawn across land masses or across the
  date-line.  This ensures sea ice and fluxes only appear over ocean cells.

- ELM regular-grid variables (QSNWCPICE, SNO_T_davg) are rendered directly with
  pcolormesh on their native 2-D lat/lon grid — no interpolation.

Usage
-----
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_login.sh
    python plot_maps.py

Edit config.py to change run directories, map extent, or output location.
"""

import os
import re
import sys

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import cartopy.crs as ccrs
import cartopy.feature as cfeature

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)

import config
import utils


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MONTH_NAMES = [
    'January', 'February', 'March',     'April',   'May',      'June',
    'July',    'August',   'September', 'October', 'November', 'December',
]


# ---------------------------------------------------------------------------
# Variable definitions
# ---------------------------------------------------------------------------
# Each entry fully describes one variable to map.  The 'source' key routes
# loading to the correct reader function and determines the rendering path.
# 'vmin'/'vmax' override auto-scaling when set (e.g. iceAreaCell → [0, 1]).

MAP_VARIABLES = [
    dict(
        name         = 'iceAreaCell',
        long_name    = 'Monthly Sea-Ice Area Fraction',
        units        = 'fraction',
        cmap         = 'Blues',
        source       = 'mpassi',
        file_pattern = '*.mpassi.hist.am.timeSeriesStatsMonthly.*.nc',
        var_cands    = ['timeMonthly_avg_iceAreaCell',
                        'timeMonthly_avg_iceArea', 'iceAreaCell'],
        depth_dim    = None,
        vmin         = 0.0,
        vmax         = 1.0,
    ),
    dict(
        name         = 'iceRunoffFlux',
        long_name    = 'Monthly Ice Runoff Flux',
        units        = 'kg m⁻² s⁻¹',
        cmap         = 'YlOrRd',
        source       = 'mpaso',
        file_pattern = '*.mpaso.hist.am.timeSeriesStatsMonthly.*.nc',
        var_cands    = ['timeMonthly_avg_iceRunoffFlux', 'iceRunoffFlux'],
        depth_dim    = None,
        vmin         = None,
        vmax         = None,
    ),
    dict(
        name         = 'riverRunoffFlux',
        long_name    = 'Monthly River Runoff Flux',
        units        = 'kg m⁻² s⁻¹',
        cmap         = 'YlGnBu',
        source       = 'mpaso',
        file_pattern = '*.mpaso.hist.am.timeSeriesStatsMonthly.*.nc',
        var_cands    = ['timeMonthly_avg_riverRunoffFlux', 'riverRunoffFlux'],
        depth_dim    = None,
        vmin         = None,
        vmax         = None,
    ),
    dict(
        name         = 'QSNWCPICE',
        long_name    = 'Snow-Cap Ice Flux (QSNWCPICE)',
        units        = 'mm s⁻¹',
        cmap         = 'viridis',
        source       = 'elm_h0',
        stream       = 'h0',
        var_name     = 'QSNWCPICE',
        depth_dim    = None,
        vmin         = None,
        vmax         = None,
    ),
    dict(
        name         = 'SNO_T_davg',
        long_name    = 'Snow Temperature — depth-averaged (SNO_T)',
        units        = 'K',
        cmap         = 'RdYlBu_r',
        source       = 'elm_h2',
        stream       = 'h2',
        var_name     = 'SNO_T',
        depth_dim    = 'levsno',   # average over all snow layers
        vmin         = None,
        vmax         = None,
    ),
]


# ---------------------------------------------------------------------------
# MPAS data loading
# ---------------------------------------------------------------------------

# _mpas_coords is no longer used for data files — mesh coords come from
# restart files via utils.load_mpassi_mesh / utils.load_mpaso_mesh.


def _month_from_mpas_path(fpath, ds=None):
    """
    Extract the month (1-12) from an MPAS analysis-member filename or, as a
    fallback, from the xtime variable inside the file.

    Handles:
      *.TYPE.YYYY-MM-DD.nc    → month = MM
      *.TYPE.YYYY.MM.nc       → month = MM
    """
    base = os.path.basename(fpath)
    # YYYY-MM-DD.nc
    m = re.search(r'\.(\d{4})-(\d{2})-\d{2}\.nc$', base)
    if m:
        return int(m.group(2))
    # YYYY.MM.nc
    m = re.search(r'\.(\d{4})\.(\d{2})\.nc$', base)
    if m:
        return int(m.group(2))
    # Fallback: read xtime from the open dataset
    if ds is not None and 'xtime' in ds:
        try:
            raw = ds['xtime'].values
            s = (''.join(raw[0].astype(str))
                 if raw.ndim == 2
                 else str(raw[0]))
            return int(s.strip()[5:7])
        except Exception:
            pass
    print(f"  Warning: cannot determine month for {base}")
    return None


def load_mpas_monthly(run_dir, file_pattern, var_cands, mesh_lat, mesh_lon):
    """
    Load a single scalar variable from MPAS monthly analysis-member files.

    mesh_lat, mesh_lon : 1-D arrays of cell coordinates in degrees, loaded
                         from the corresponding component restart file (these
                         variables are absent from the analysis-member history
                         files to save disk space).

    Returns
    -------
    dict : month_int (1–12) → (lat_deg 1-D, lon_deg 1-D, data 1-D)
    Works for multi-year runs: all years are loaded; where a month appears in
    multiple years the last year's value overwrites earlier ones (single-year
    climatology not yet implemented — values are taken as-is).
    """
    files = utils.find_files(run_dir, file_pattern)
    if not files:
        return {}

    result = {}
    printed_var = False

    for fpath in files:
        try:
            ds = xr.open_dataset(fpath, decode_times=False, mask_and_scale=True)
        except Exception as exc:
            print(f"  Warning: cannot open {os.path.basename(fpath)}: {exc}")
            continue

        month = _month_from_mpas_path(fpath, ds)
        if month is None:
            ds.close()
            continue

        # Find the requested variable
        var_name = None
        for cand in var_cands:
            if cand in ds:
                var_name = cand
                break
        if var_name is None:
            if not printed_var:
                avail = list(ds.data_vars)
                print(f"  Warning: none of {var_cands} found. "
                      f"Available (first 30): {avail[:30]}")
                printed_var = True
            ds.close()
            continue

        data = ds[var_name].values
        # Squeeze Time/time dimension (one record per file)
        for dim in ('Time', 'time'):
            if dim in ds[var_name].dims:
                data = data[0]
                break
        data = data.ravel().astype(float)
        data = np.where(np.abs(data) < 1e30, data, np.nan)

        result[month] = (mesh_lat, mesh_lon, data)
        ds.close()

    n = len(result)
    if n:
        print(f"  Loaded {n} month(s) from {os.path.basename(file_pattern)}")
    return result


# ---------------------------------------------------------------------------
# ELM data loading
# ---------------------------------------------------------------------------

def _elm_coords(ds):
    """Return (lat_deg, lon_deg) from an ELM xarray Dataset."""
    for v in ['lat', 'LATIXY', 'latitude']:
        if v in ds:
            lat = ds[v].values.ravel()
            break
    else:
        raise KeyError(f"No latitude variable found in ELM file. "
                       f"Available: {list(ds.coords) + list(ds.data_vars)}")
    for v in ['lon', 'LONGXY', 'longitude']:
        if v in ds:
            lon = ds[v].values.ravel()
            break
    else:
        raise KeyError("No longitude variable found in ELM file.")
    lon = utils.normalize_lon(lon)
    return lat, lon


def _elm_decode_months(ds, tdim):
    """
    Return {time_index: month_int} for every time step in *ds* by decoding
    the ELM time coordinate (units / calendar attributes).

    ELM monthly history files often contain all 12 months in a single file
    (time values = days since reference, noleap calendar).  cftime is used
    to decode the values; returns an empty dict if decoding fails.
    """
    try:
        import cftime
        time_var = ds[tdim]
        units    = time_var.attrs.get('units', '')
        calendar = time_var.attrs.get('calendar', 'standard')
        dates    = cftime.num2date(time_var.values, units, calendar)
        return {i: int(d.month) for i, d in enumerate(dates)}
    except Exception as exc:
        print(f"  Warning: could not decode ELM time axis ({exc})")
        return {}


def load_elm_monthly(run_dir, stream, var_name, depth_dim=None):
    """
    Load a monthly ELM variable from h0 or h2 history files.

    Handles both the common case where each file holds a single month and
    the case where a single file contains all 12 months (time dim > 1).
    The month for each record is determined by decoding the time coordinate
    via cftime; the filename-derived month is used as a fallback for
    single-record files.

    depth_dim : if not None, the named dimension is averaged over with
                nanmean (used for SNO_T levsno depth-averaging).
                Works for any number of snow layers (standard 5-layer or
                extended snowpack).

    Returns
    -------
    dict : month_int (1–12) → (lat_1d, lon_1d, data_2d)
           lat_1d and lon_1d are 1-D coordinate axes; data_2d is shape
           (nlat, nlon).  This preserves the native grid structure for
           direct pcolormesh rendering without interpolation.
    """
    by_month = utils.get_elm_files_by_month(run_dir, stream)
    if not by_month:
        return {}

    result = {}
    printed_var = False

    for (year, fname_month), fpath in sorted(by_month.items()):
        try:
            ds = xr.open_dataset(fpath, decode_times=False, mask_and_scale=True)
        except Exception as exc:
            print(f"  Warning: cannot open {os.path.basename(fpath)}: {exc}")
            continue

        if var_name not in ds:
            if not printed_var:
                avail = list(ds.data_vars)
                print(f"  Warning: '{var_name}' not found in elm.{stream} files. "
                      f"Available (first 30): {avail[:30]}")
                printed_var = True
            ds.close()
            continue

        try:
            lat_raw, lon_raw = _elm_coords(ds)
        except KeyError as exc:
            print(f"  Warning: {exc}")
            ds.close()
            continue

        da = ds[var_name]

        # Identify the time dimension (if any)
        tdim = None
        for d in ('time', 'Time'):
            if d in da.dims:
                tdim = d
                break

        # Depth-average before iterating over time (e.g. SNO_T over levsno)
        if depth_dim is not None:
            if depth_dim in da.dims:
                da = da.mean(dim=depth_dim, skipna=True)
            else:
                for candidate in ['levsno', 'levtot', 'levsnl', 'snowlayer']:
                    if candidate in da.dims:
                        if not printed_var:
                            print(f"  Note: '{depth_dim}' not found; "
                                  f"averaging over '{candidate}' instead.")
                        da = da.mean(dim=candidate, skipna=True)
                        break
                else:
                    if not printed_var:
                        print(f"  Warning: depth dim '{depth_dim}' not found; "
                              f"using data as-is (dims={list(da.dims)}).")

        # Map each time index to its calendar month.
        # When a file holds multiple months, decode from the time coordinate.
        if tdim is not None and da.sizes[tdim] > 1:
            ti_to_month = _elm_decode_months(ds, tdim)
        elif tdim is not None:
            ti_to_month = {0: fname_month}
        else:
            ti_to_month = {None: fname_month}   # no time dimension at all

        for ti, month in sorted(ti_to_month.items()):
            if month is None:
                continue

            # Extract this time step (or use the whole array if no time dim)
            da_t = da.isel({tdim: ti}) if (tdim is not None and ti is not None) else da

            data = da_t.values.astype(float)
            data = np.where(np.abs(data) < 1e30, data, np.nan)

            # Ensure data is 2-D (nlat, nlon).  If the ELM file uses a flat
            # gridcell dimension, reshape using the lat/lon axis lengths.
            if data.ndim == 1 and lat_raw.size != data.size:
                data = data.reshape(len(lat_raw), len(lon_raw))
            elif data.ndim == 1:
                # Unstructured — cannot preserve 2D; reshape if possible
                nlat = len(np.unique(lat_raw))
                nlon = data.size // nlat if nlat > 0 else 0
                if nlat * nlon == data.size:
                    data = data.reshape(nlat, nlon)
                    lat_raw = np.unique(lat_raw)
                    lon_raw = np.unique(lon_raw)

            result[month] = (lat_raw.copy(), lon_raw.copy(), data)

        ds.close()

    return result


# ---------------------------------------------------------------------------
# MPAS triangulation helper
# ---------------------------------------------------------------------------

def _build_masked_triangulation(lat, lon):
    """
    Build a matplotlib Triangulation from MPAS cell centers (lat, lon in
    degrees) with a mask that removes:
      - Triangles spanning the date-line (huge edge in longitude)
      - Triangles spanning land masses (edges much larger than typical cell
        spacing)

    The edge-length threshold is set to 2× the median edge length of the
    Delaunay triangulation, which adapts automatically to any mesh resolution.

    Returns a matplotlib.tri.Triangulation object (reusable across all months
    for the same mesh).
    """
    tri = mtri.Triangulation(lon, lat)
    triangles = tri.triangles  # (N_tri, 3) indices

    # Compute edge lengths (in degrees) for each triangle
    x, y = lon, lat
    x0 = x[triangles[:, 0]]; x1 = x[triangles[:, 1]]; x2 = x[triangles[:, 2]]
    y0 = y[triangles[:, 0]]; y1 = y[triangles[:, 1]]; y2 = y[triangles[:, 2]]

    # Compute squared edge lengths (sufficient for comparison)
    dx01 = x1 - x0; dy01 = y1 - y0
    dx12 = x2 - x1; dy12 = y2 - y1
    dx20 = x0 - x2; dy20 = y0 - y2

    edge2_01 = dx01**2 + dy01**2
    edge2_12 = dx12**2 + dy12**2
    edge2_20 = dx20**2 + dy20**2

    max_edge2 = np.maximum(np.maximum(edge2_01, edge2_12), edge2_20)

    # Threshold: 2× median edge length (squared for comparison)
    all_edge2 = np.concatenate([edge2_01, edge2_12, edge2_20])
    median_edge = np.sqrt(np.median(all_edge2))
    threshold = (2.0 * median_edge) ** 2

    # Mask triangles with any edge exceeding threshold
    mask = max_edge2 > threshold
    tri.set_mask(mask)

    n_masked = mask.sum()
    n_total  = len(triangles)
    print(f"    Triangulation: {len(lat)} points, {n_total} triangles, "
          f"{n_masked} masked ({100*n_masked/n_total:.1f}%)")
    return tri


# ---------------------------------------------------------------------------
# Colour scale helpers
# ---------------------------------------------------------------------------

def _compute_vrange(all_data_arrays, fixed_vmin=None, fixed_vmax=None):
    """
    Return (vmin, vmax) for colour scaling.

    If fixed_vmin/vmax are provided (not None), those values are used directly.
    Otherwise, computes the 2nd/98th percentile from all non-NaN values in the
    provided arrays (which may be 1-D or 2-D).

    Fallback for very sparse fields (where p2 == p98, e.g. >98% of values are
    zero): uses the actual min/max of non-zero values, with a floor at 0 for
    non-negative data.
    """
    if fixed_vmin is not None and fixed_vmax is not None:
        return float(fixed_vmin), float(fixed_vmax)

    vals = []
    for arr in all_data_arrays:
        if arr is not None:
            flat = np.asarray(arr).ravel()
            finite = flat[np.isfinite(flat)]
            if len(finite) > 0:
                vals.append(finite)
    if not vals:
        return 0.0, 1.0
    all_vals = np.concatenate(vals)
    vmin = float(np.nanpercentile(all_vals, 2))  if fixed_vmin is None else float(fixed_vmin)
    vmax = float(np.nanpercentile(all_vals, 98)) if fixed_vmax is None else float(fixed_vmax)
    if vmin == vmax:
        # Sparse field: percentiles are identical (often both zero).
        # Use the actual data range of non-zero values instead.
        nonzero = all_vals[all_vals != 0.0]
        if len(nonzero) > 0:
            data_min = float(nonzero.min())
            data_max = float(nonzero.max())
            # For non-negative data, floor at 0
            vmin = 0.0 if data_min >= 0 else data_min
            vmax = data_max
        else:
            # Truly all zeros
            vmin, vmax = 0.0, 1.0
        # Small guard against vmin == vmax after nonzero range
        if vmin == vmax:
            vmin, vmax = vmin - 0.5 * abs(vmin + 1e-10), vmax + 0.5 * abs(vmax + 1e-10)
    return vmin, vmax


# ---------------------------------------------------------------------------
# Per-variable map generation
# ---------------------------------------------------------------------------

def _render_panel(ax, lat_d, lon_d, data_d, source_type, tri_cache,
                  cache_key, var_name, proj, cmap, vmin, vmax):
    """
    Render data onto a single map axes panel.  Returns the mappable or None.
    """
    if source_type == 'unstructured':
        # --- MPAS: tripcolor with masked triangulation ---
        tri_key = (cache_key, var_name, id(lat_d))
        if tri_key not in tri_cache:
            tri_cache[tri_key] = _build_masked_triangulation(lat_d, lon_d)
        tri = tri_cache[tri_key]

        plot_data = data_d.copy()
        nan_mask = ~np.isfinite(plot_data)
        plot_data[nan_mask] = 0.0

        tri_verts = tri.triangles
        all_nan = (nan_mask[tri_verts[:, 0]]
                   & nan_mask[tri_verts[:, 1]]
                   & nan_mask[tri_verts[:, 2]])
        combined_mask = tri.mask | all_nan if tri.mask is not None else all_nan

        tri_local = mtri.Triangulation(tri.x, tri.y, tri.triangles)
        tri_local.set_mask(combined_mask)

        im = ax.tripcolor(tri_local, plot_data,
                          transform=proj,
                          cmap=cmap,
                          vmin=vmin, vmax=vmax,
                          shading='flat',
                          zorder=1)
        return im

    else:
        # --- ELM: direct pcolormesh on native 2D grid ---
        if data_d.ndim == 2 and len(lat_d) == data_d.shape[0] \
                and len(lon_d) == data_d.shape[1]:
            im = ax.pcolormesh(lon_d, lat_d, data_d,
                               transform=proj,
                               cmap=cmap,
                               vmin=vmin, vmax=vmax,
                               shading='auto',
                               zorder=1)
        else:
            flat = data_d.ravel()
            lat_f = lat_d if lat_d.size == flat.size else \
                np.meshgrid(lat_d, lon_d, indexing='ij')[0].ravel()
            lon_f = lon_d if lon_d.size == flat.size else \
                np.meshgrid(lat_d, lon_d, indexing='ij')[1].ravel()
            im = ax.scatter(lon_f, lat_f, c=flat, s=0.5,
                            transform=proj, cmap=cmap,
                            vmin=vmin, vmax=vmax, zorder=1)
        return im


def plot_variable(var_name, long_name, units, cmap,
                  fix_data, ctrl_data, source_type,
                  lat_min, lat_max, lon_min, lon_max,
                  output_dir, tri_cache,
                  fixed_vmin=None, fixed_vmax=None,
                  show_diff=False, diff_cmap='RdBu_r'):
    """
    Generate one PNG per month for the given variable.

    Layout:
      - Standard (show_diff=False): 2 rows × 1 col
            Row 1: Control (no fix)
            Row 2: Fix run
      - With difference (show_diff=True): 3 rows × 1 col
            Row 1: Control (no fix)
            Row 2: Fix run
            Row 3: Difference (fix − control)

    A shared colourbar is placed below the map panels in its own dedicated
    space.  The difference panel gets its own symmetric colourbar.

    Parameters
    ----------
    fix_data, ctrl_data : dicts  month_int → tuple
        For source_type='unstructured' (MPAS):
            tuple = (lat_1d, lon_1d, data_1d)
        For source_type='gridded' (ELM):
            tuple = (lat_1d, lon_1d, data_2d)
    source_type : 'unstructured' or 'gridded'
    tri_cache : dict — built lazily for MPAS triangulations
    fixed_vmin, fixed_vmax : optional hard colour bounds
    show_diff : if True, add a third panel showing fix − control
    diff_cmap : colormap for the difference panel
    """
    out_dir = os.path.join(output_dir, 'maps', var_name)
    os.makedirs(out_dir, exist_ok=True)

    all_months = sorted(set(list(fix_data.keys()) + list(ctrl_data.keys())))
    if not all_months:
        print(f"  No data found for {var_name}; skipping.")
        return

    # Collect all data arrays for colour scaling
    all_data = []
    for month_data in [fix_data, ctrl_data]:
        for month in month_data:
            _, _, d = month_data[month]
            all_data.append(d)

    vmin, vmax = _compute_vrange(all_data, fixed_vmin, fixed_vmax)

    # Determine difference colour range (symmetric about zero)
    diff_vmax_abs = None
    if show_diff:
        diff_vals = []
        for month in all_months:
            if month in fix_data and month in ctrl_data:
                _, _, d_fix = fix_data[month]
                _, _, d_ctrl = ctrl_data[month]
                if d_fix.shape == d_ctrl.shape:
                    diff = d_fix - d_ctrl
                    finite = diff[np.isfinite(diff)].ravel()
                    if len(finite) > 0:
                        diff_vals.append(finite)
        if diff_vals:
            all_diff = np.concatenate(diff_vals)
            diff_vmax_abs = float(np.nanpercentile(np.abs(all_diff), 98))
            if diff_vmax_abs == 0:
                diff_vmax_abs = 1.0

    proj   = ccrs.PlateCarree()
    extent = [lon_min, lon_max, lat_min, lat_max]

    nrows = 3 if show_diff else 2

    for month in all_months:
        mname = _MONTH_NAMES[month - 1]

        # Use gridspec for explicit control of colourbar space
        fig = plt.figure(figsize=(12, 4.0 * nrows + 0.8))
        if show_diff:
            # 3 map rows + 2 colourbar rows (interleaved):
            # Row 0: Control, Row 1: Fix, Row 2: main cbar,
            # Row 3: Difference, Row 4: diff cbar
            gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 0.06, 1, 0.06],
                                  hspace=0.40,
                                  top=0.95, bottom=0.03)
        else:
            # 2 map rows + 1 colourbar row
            gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.06],
                                  hspace=0.30)

        # --- Panel definitions: (row_idx, label, data_dict, cache_key) ---
        panels = [
            (0, config.CTRL_LABEL, ctrl_data, 'ctrl'),
            (1, config.FIX_LABEL,  fix_data,  'fix'),
        ]

        mappable = None
        for row_idx, sim_label, sim_data, cache_key in panels:
            ax = fig.add_subplot(gs[row_idx, 0], projection=proj)
            ax.set_extent(extent, crs=proj)
            ax.add_feature(cfeature.LAND,      facecolor='#d8d8d8', zorder=0)
            ax.add_feature(cfeature.OCEAN,     facecolor='#e8f4f8', zorder=0)
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5,       zorder=2)
            ax.add_feature(cfeature.BORDERS,   linewidth=0.3,       zorder=2)
            gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5,
                              x_inline=False, y_inline=False)
            gl.top_labels    = False
            gl.right_labels  = False

            if month not in sim_data:
                ax.text(0.5, 0.5, 'no data',
                        transform=ax.transAxes,
                        ha='center', va='center',
                        fontsize=11, color='grey')
            else:
                lat_d, lon_d, data_d = sim_data[month]
                im = _render_panel(ax, lat_d, lon_d, data_d, source_type,
                                   tri_cache, cache_key, var_name, proj,
                                   cmap, vmin, vmax)
                if im is not None and mappable is None:
                    mappable = im

            ax.set_title(f'{sim_label}  —  {mname}', fontsize=10)

        # --- Colourbar for the two main panels ---
        cbar_row = 2 if show_diff else 2
        cax = fig.add_subplot(gs[cbar_row, 0])
        if mappable is not None:
            cbar = fig.colorbar(mappable, cax=cax, orientation='horizontal')
            cbar.set_label(f'{long_name}  ({units})', fontsize=9)
        else:
            cax.set_visible(False)

        # --- Difference panel (row 3, if show_diff) ---
        diff_mappable = None
        if show_diff:
            ax_diff = fig.add_subplot(gs[3, 0], projection=proj)
            ax_diff.set_extent(extent, crs=proj)
            ax_diff.add_feature(cfeature.LAND,      facecolor='#d8d8d8', zorder=0)
            ax_diff.add_feature(cfeature.OCEAN,     facecolor='#e8f4f8', zorder=0)
            ax_diff.add_feature(cfeature.COASTLINE, linewidth=0.5,       zorder=2)
            ax_diff.add_feature(cfeature.BORDERS,   linewidth=0.3,       zorder=2)
            gl = ax_diff.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5,
                                   x_inline=False, y_inline=False)
            gl.top_labels    = False
            gl.right_labels  = False

            if month in fix_data and month in ctrl_data:
                lat_f, lon_f, d_fix = fix_data[month]
                lat_c, lon_c, d_ctrl = ctrl_data[month]

                if d_fix.shape == d_ctrl.shape:
                    diff_data = d_fix - d_ctrl
                    dvm = diff_vmax_abs if diff_vmax_abs else 1.0
                    im_diff = _render_panel(
                        ax_diff, lat_f, lon_f, diff_data, source_type,
                        tri_cache, 'fix', var_name, proj,
                        diff_cmap, -dvm, dvm)
                    diff_mappable = im_diff
                else:
                    ax_diff.text(0.5, 0.5, 'shape mismatch',
                                transform=ax_diff.transAxes,
                                ha='center', va='center',
                                fontsize=11, color='grey')
            else:
                ax_diff.text(0.5, 0.5, 'no data for diff',
                            transform=ax_diff.transAxes,
                            ha='center', va='center',
                            fontsize=11, color='grey')

            ax_diff.set_title(
                f'Difference ({config.FIX_LABEL} − {config.CTRL_LABEL})  —  {mname}',
                fontsize=10)

            # Difference colourbar
            cax_diff = fig.add_subplot(gs[4, 0])
            if diff_mappable is not None:
                cbar_diff = fig.colorbar(diff_mappable, cax=cax_diff,
                                         orientation='horizontal')
                cbar_diff.set_label(f'Δ {long_name}  ({units})', fontsize=9)
            else:
                cax_diff.set_visible(False)

        fig.suptitle(f'{long_name}  |  {mname}', fontsize=12, y=0.99)

        outpath = os.path.join(out_dir, f'month_{month:02d}_{var_name}.png')
        fig.savefig(outpath, dpi=config.DPI, bbox_inches='tight')
        plt.close(fig)

    print(f"  Saved {len(all_months)} map(s) for '{var_name}' → {out_dir}/")


# ---------------------------------------------------------------------------
# Loader dispatcher
# ---------------------------------------------------------------------------

def load_variable(run_dir, var_def, mesh_cache):
    """
    Route loading to the correct reader based on var_def['source'].

    mesh_cache : dict keyed by (run_dir, component) → (lat, lon) tuples for
                 each MPAS component; populated lazily on first access.

    Returns
    -------
    (data_dict, source_type)
        data_dict   : month_int → tuple (see load_mpas_monthly / load_elm_monthly)
        source_type : 'unstructured' (MPAS) or 'gridded' (ELM)
    """
    src = var_def['source']
    if src in ('mpassi', 'mpaso'):
        # Load mesh coords from restart file (once per run per component)
        cache_key = (run_dir, src)
        if cache_key not in mesh_cache:
            try:
                if src == 'mpassi':
                    lat, lon, _ = utils.load_mpassi_mesh(run_dir)
                else:
                    lat, lon, _ = utils.load_mpaso_mesh(run_dir)
                mesh_cache[cache_key] = (lat, lon)
            except RuntimeError as exc:
                print(f"  Error loading mesh: {exc}")
                return {}, 'unstructured'
        mesh_lat, mesh_lon = mesh_cache[cache_key]
        data = load_mpas_monthly(
            run_dir,
            var_def['file_pattern'],
            var_def['var_cands'],
            mesh_lat,
            mesh_lon,
        )
        return data, 'unstructured'
    elif src in ('elm_h0', 'elm_h2'):
        data = load_elm_monthly(
            run_dir,
            var_def['stream'],
            var_def['var_name'],
            depth_dim=var_def.get('depth_dim'),
        )
        return data, 'gridded'
    print(f"  Warning: unknown source '{src}' for variable '{var_def['name']}'")
    return {}, 'unstructured'


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Monthly Map Plots")
    print("=" * 60)
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)

    lat_min, lat_max, lon_min, lon_max = utils.resolve_map_extent(
        config.MAP_LAT_MIN, config.MAP_LAT_MAX,
        config.MAP_LON_MIN, config.MAP_LON_MAX,
        config.FIX_RUN_DIR,
    )

    # Mesh coordinate cache: populated lazily the first time each
    # (run_dir, component) pair is encountered, then reused across variables.
    mesh_cache = {}

    # Triangulation cache: populated lazily per (run_label, var_name, mesh_id);
    # reused across months for the same mesh.
    tri_cache = {}

    # Filter variables if user specified a subset in config
    vars_to_plot = getattr(config, 'MAP_VARIABLES_TO_PLOT', None)
    if vars_to_plot is not None:
        active_vars = [v for v in MAP_VARIABLES if v['name'] in vars_to_plot]
        skipped = [v['name'] for v in MAP_VARIABLES if v['name'] not in vars_to_plot]
        if skipped:
            print(f"  Skipping (per MAP_VARIABLES_TO_PLOT): {', '.join(skipped)}")
    else:
        active_vars = MAP_VARIABLES

    for var_def in active_vars:
        vname = var_def['name']
        print(f"\n{'─' * 50}")
        print(f"Variable: {vname}")

        print(f"  Loading fix run…")
        fix_data, source_type = load_variable(config.FIX_RUN_DIR, var_def, mesh_cache)
        print(f"  Loading ctrl run…")
        ctrl_data, _ = load_variable(config.CTRL_RUN_DIR, var_def, mesh_cache)

        if not fix_data and not ctrl_data:
            print(f"  No data found for {vname}; skipping.")
            continue

        # Determine colour range: user override > var_def > auto
        vrange_override = getattr(config, 'MAP_VRANGE_OVERRIDES', {}).get(vname)
        if vrange_override is not None:
            fv_min, fv_max = vrange_override
        else:
            fv_min = var_def.get('vmin')
            fv_max = var_def.get('vmax')

        plot_variable(
            var_name    = vname,
            long_name   = var_def['long_name'],
            units       = var_def['units'],
            cmap        = var_def['cmap'],
            fix_data    = fix_data,
            ctrl_data   = ctrl_data,
            source_type = source_type,
            lat_min=lat_min, lat_max=lat_max,
            lon_min=lon_min, lon_max=lon_max,
            output_dir  = config.OUTPUT_DIR,
            tri_cache   = tri_cache,
            fixed_vmin  = fv_min,
            fixed_vmax  = fv_max,
            show_diff   = (vname == 'SNO_T_davg'),
        )

    print("\nDone.\n")


if __name__ == '__main__':
    main()
