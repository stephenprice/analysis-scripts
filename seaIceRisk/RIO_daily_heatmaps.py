#!/usr/bin/env python3
"""
RIO_daily_heatmaps.py
=====================
Standalone script version of RIO_daily_heatmaps.ipynb.

Produces heatmap figures (day-of-year × year) of POLARIS sea ice shipping
risk for one or more vessel classes and Arctic route sections, from daily
MPAS-SI ensemble-statistics NetCDF files.

Usage
-----
    python RIO_daily_heatmaps.py

Edit the USER CONFIGURATION section below to change paths, vessel classes,
route, years, plot mode, etc.
"""

import os
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')          # non-interactive backend; figures are saved to disk
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from scipy.spatial import KDTree
from datetime import date


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# --- Paths ---
data_path  = os.path.expanduser('~/projects/CINS/ColdHarbor/test-data/')
mesh_file  = os.path.join(data_path, 'mpaso-IcoswISC30E3r5-restart.nc')
route_file = os.path.expanduser(
    '~/projects/CINS/ColdHarbor/pythonScripts/arctic_route_sections.nc')
output_dir = os.path.expanduser('~/projects/CINS/ColdHarbor/pythonScripts/')

# --- Vessel types to process ---
# Available: 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',
#            'IA', 'IAsuper', 'IB', 'IC', 'NIS'
vessel_classes = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',
                  'IA', 'IAsuper', 'IB', 'IC', 'NIS']

# --- Route sections to concatenate, in order ---
# Available sections: A through J (as defined in arctic_route_sections.nc)
route_sections = ['E', 'J', 'C', 'D']
# route_sections = ['H', 'I', 'B', 'A']

# --- Years to include ---
year_start = 2000
year_end   = 2019
years = [str(y) for y in range(year_start, year_end + 1)]

# --- Plot mode ---
# 'transit_time' : segment-integrated transit time in days (NaN = impassable)
# 'RIO'          : route RIO summary value (method set by rio_summary_method)
# 'passability'  : binary green (passable) / red (impassable)
plot_mode = 'transit_time'
# plot_mode = 'RIO'
# plot_mode = 'passability'

# --- Output ---
save_figs = True    # set False to suppress file saving

# --- RIO heatmap summary method ---
# Controls what value is shown in RIO heatmaps. Passability always uses min.
# 'min'  : route-minimum RIO (worst bottleneck along the route)
# 'mean' : route-mean RIO (average ice severity experienced along the transit)
rio_summary_method = 'min'

# --- x-axis centering ---
# 0-based day-of-year to place at the centre of the heatmap x-axis.
# 258 ~ September 16, near the annual Arctic sea-ice minimum.
plot_center_doy = 258

# --- Spatial interpolation ---
# Number of nearest MPAS mesh cells to average RIO over at each route waypoint
k_neighbors = 4


# =============================================================================
# VESSEL SPEED / IMPASSABILITY CONFIGURATION
# =============================================================================

VESSEL_CONFIG = {
    'PC1':     {'normal_kt': 19.0,  'restricted_kt': 11.0, 'impass_thresh': -10},
    'PC2':     {'normal_kt': 16.5,  'restricted_kt':  8.0, 'impass_thresh': -10},
    'PC3':     {'normal_kt': 16.0,  'restricted_kt':  5.0, 'impass_thresh': -10},
    'PC4':     {'normal_kt': 15.0,  'restricted_kt':  5.0, 'impass_thresh': -10},
    'PC5':     {'normal_kt': 14.0,  'restricted_kt':  5.0, 'impass_thresh': -10},
    'PC6':     {'normal_kt': 14.0,  'restricted_kt':  3.0, 'impass_thresh': -10},
    'PC7':     {'normal_kt': 14.0,  'restricted_kt':  3.0, 'impass_thresh': -10},
    'IAsuper': {'normal_kt': 14.0,  'restricted_kt': None, 'impass_thresh':   0},
    'IA':      {'normal_kt': 14.0,  'restricted_kt': None, 'impass_thresh':   0},
    'IB':      {'normal_kt': 14.0,  'restricted_kt': None, 'impass_thresh':   0},
    'IC':      {'normal_kt': 14.0,  'restricted_kt': None, 'impass_thresh':   0},
    'NIS':     {'normal_kt': 14.0,  'restricted_kt': None, 'impass_thresh':   0},
}


def rio_to_speed_array(rio_arr, vessel_class):
    """
    Map an array of per-waypoint RIO values to an array of speeds (knots).

    Waypoints that are impassable (RIO < impass_thresh) receive NaN.
    Non-PC vessels have no restricted band: they sail at normal_kt or are
    impassable.

    Parameters
    ----------
    rio_arr      : ndarray, shape (n_waypoints,)
    vessel_class : str   key into VESSEL_CONFIG

    Returns
    -------
    speed : ndarray, shape (n_waypoints,)  in knots; NaN where impassable
    """
    cfg    = VESSEL_CONFIG[vessel_class]
    speed  = np.full(len(rio_arr), np.nan, dtype=float)
    thresh = cfg['impass_thresh']

    passable   = (~np.isnan(rio_arr)) & (rio_arr >= thresh)
    normal     = passable & (rio_arr >= 0)
    restricted = passable & (rio_arr < 0)

    speed[normal] = cfg['normal_kt']
    if cfg['restricted_kt'] is not None:
        speed[restricted] = cfg['restricted_kt']
    else:
        # non-PC: no restricted band — sail at normal speed when passable
        speed[restricted] = cfg['normal_kt']

    return speed


def integrate_transit_time(wp_rio, seg_dist_nm, vessel_class):
    """
    Compute segment-integrated transit time in days.

    Each waypoint is assigned a local speed based on its RIO value.
    Transit time = sum(seg_dist_i / speed_i) / 24.
    Returns NaN if any waypoint along the route is impassable.

    Parameters
    ----------
    wp_rio        : ndarray (n_waypoints,)  per-waypoint RIO (k-neighbour average)
    seg_dist_nm   : ndarray (n_waypoints,)  distance represented by each waypoint (NM)
    vessel_class  : str

    Returns
    -------
    float : transit time in days, or NaN if route is impassable
    """
    speed = rio_to_speed_array(wp_rio, vessel_class)
    if np.any(np.isnan(speed)):
        return np.nan
    return float(np.sum(seg_dist_nm / speed) / 24.0)


# =============================================================================
# LOAD MPAS MESH
# =============================================================================

print('Loading MPAS mesh...')
ds_mesh = xr.open_dataset(mesh_file)

rad2deg       = 180.0 / np.pi
lat_mesh_full = ds_mesh.latCell.values * rad2deg
lon_mesh_full = ds_mesh.lonCell.values * rad2deg

# Normalise longitudes to [-180, 180]
lon_mesh_full = np.where(lon_mesh_full > 180.0,
                         lon_mesh_full - 360.0,
                         lon_mesh_full)

# Restrict to cells at or north of 60 deg N
ind_north = np.where(lat_mesh_full >= 60.0)[0]
lat_mesh  = lat_mesh_full[ind_north]
lon_mesh  = lon_mesh_full[ind_north]

print(f'  Total mesh cells  : {len(lat_mesh_full):,}')
print(f'  Cells >= 60 deg N : {len(lat_mesh):,}')


# =============================================================================
# LOAD ROUTE SECTIONS AND BUILD KDTREE
# =============================================================================

print('Loading route sections...')
ds_route = xr.open_dataset(route_file)

trans_lat_list = []
trans_lon_list = []
seg_dist_list  = []   # differential distance per waypoint (NM)

for sec in route_sections:
    s   = sec.upper()
    dis = ds_route[f'dis_{s}'].values   # cumulative distance (NM), starts at 0
    trans_lat_list.append(ds_route[f'lat_{s}'].values)
    trans_lon_list.append(ds_route[f'lon_{s}'].values)
    # np.diff gives the inter-waypoint spacing; prepend=0 so that the first
    # waypoint of each section contributes dis[0] NM (which is 0 by convention)
    seg_dist_list.append(np.diff(dis, prepend=0.0))

trans_lat     = np.concatenate(trans_lat_list)
trans_lon     = np.concatenate(trans_lon_list)
seg_dist_nm   = np.concatenate(seg_dist_list)   # (n_waypoints,)
total_dist_nm = sum(
    float(np.max(ds_route[f'dis_{s.upper()}'].values)) for s in route_sections
)
route_label = ''.join(s.upper() for s in route_sections)

print(f'  Sections        : {route_sections}  ->  {len(trans_lat)} waypoints')
print(f'  Total distance  : {total_dist_nm:.1f} nm')
print(f'  seg_dist_nm     : min={seg_dist_nm.min():.3f}, '
      f'max={seg_dist_nm.max():.3f}, sum={seg_dist_nm.sum():.1f} nm')

print('Building KDTree...')
tree = KDTree(list(zip(lon_mesh, lat_mesh)))
_, inds = tree.query(list(zip(trans_lon, trans_lat)), k=k_neighbors)
print(f'  KDTree ready.  Neighbor index array shape: {inds.shape}')


# =============================================================================
# DATA EXTRACTION LOOP
# =============================================================================

def _doy_index(year, month, day):
    """Return 0-based day-of-year index (0 = Jan 1)."""
    return (date(year, month, day) - date(year, 1, 1)).days


def _find_files(vessel_class, years, data_path):
    """Return sorted list of (year, month, filepath) tuples for available monthly files."""
    prefix  = (f'v3.LR.historical_{vessel_class}_EnsStats'
               f'.mpassi.hist.am.timeSeriesStatsDaily.')
    records = []
    for yr in years:
        for mo in range(1, 13):
            fname = os.path.join(data_path, f'{prefix}{yr}-{mo:02d}-01.nc')
            if os.path.exists(fname):
                records.append((int(yr), mo, fname))
    return records


# results[vessel_class] dict holds:
#   'years'                    : list of int year values actually found
#   'min_rio_med/5th/95th'     : route-minimum RIO (used for passability; always min)
#   'summary_rio_med/5th/95th' : route RIO per rio_summary_method (used for RIO heatmap)
#   'tt_med/5th/95th'          : segment-integrated transit time (days)
results = {}

for vessel_class in vessel_classes:
    print(f'\n===== {vessel_class} =====')
    file_records = _find_files(vessel_class, years, data_path)

    if not file_records:
        print(f'  WARNING: no files found for {vessel_class} -- skipping.')
        continue

    avail_years  = sorted(set(yr for yr, _, __ in file_records))
    n_years      = len(avail_years)
    year_idx_map = {yr: i for i, yr in enumerate(avail_years)}

    shape            = (n_years, 366)
    min_rio_med      = np.full(shape, np.nan)
    min_rio_5th      = np.full(shape, np.nan)
    min_rio_95th     = np.full(shape, np.nan)
    summary_rio_med  = np.full(shape, np.nan)
    summary_rio_5th  = np.full(shape, np.nan)
    summary_rio_95th = np.full(shape, np.nan)
    tt_med           = np.full(shape, np.nan)
    tt_5th           = np.full(shape, np.nan)
    tt_95th          = np.full(shape, np.nan)

    for yr, mo, fpath in file_records:
        print(f'  {os.path.basename(fpath)}')
        ds = xr.open_dataset(fpath)

        # RIO files are already subsetted to Arctic cells; inds index directly into this space
        rio_med_arr  = ds['timeDaily_avg_RIO_ensembleMedian'].values
        rio_5th_arr  = ds['timeDaily_avg_RIO_ensemble5th'].values
        rio_95th_arr = ds['timeDaily_avg_RIO_ensemble95th'].values
        ds.close()

        n_days = rio_med_arr.shape[0]
        yi     = year_idx_map[yr]

        for day_idx in range(n_days):
            doy = _doy_index(yr, mo, day_idx + 1)   # 0-based day-of-year

            # Average k neighbours at each route waypoint: shape (n_waypoints,)
            wp_rio_med  = np.mean(rio_med_arr [day_idx, inds], axis=1)
            wp_rio_5th  = np.mean(rio_5th_arr [day_idx, inds], axis=1)
            wp_rio_95th = np.mean(rio_95th_arr[day_idx, inds], axis=1)

            # Route-minimum RIO (always stored; used for passability)
            min_rio_med [yi, doy] = float(np.nanmin(wp_rio_med))
            min_rio_5th [yi, doy] = float(np.nanmin(wp_rio_5th))
            min_rio_95th[yi, doy] = float(np.nanmin(wp_rio_95th))

            # Summary RIO (method-controlled; used for RIO heatmap)
            if rio_summary_method == 'mean':
                summary_rio_med [yi, doy] = float(np.nanmean(wp_rio_med))
                summary_rio_5th [yi, doy] = float(np.nanmean(wp_rio_5th))
                summary_rio_95th[yi, doy] = float(np.nanmean(wp_rio_95th))
            else:  # 'min'
                summary_rio_med [yi, doy] = float(np.nanmin(wp_rio_med))
                summary_rio_5th [yi, doy] = float(np.nanmin(wp_rio_5th))
                summary_rio_95th[yi, doy] = float(np.nanmin(wp_rio_95th))

            # Segment-integrated transit time
            tt_med [yi, doy] = integrate_transit_time(wp_rio_med,  seg_dist_nm, vessel_class)
            tt_5th [yi, doy] = integrate_transit_time(wp_rio_5th,  seg_dist_nm, vessel_class)
            tt_95th[yi, doy] = integrate_transit_time(wp_rio_95th, seg_dist_nm, vessel_class)

    results[vessel_class] = {
        'years':             avail_years,
        'min_rio_med':       min_rio_med,
        'min_rio_5th':       min_rio_5th,
        'min_rio_95th':      min_rio_95th,
        'summary_rio_med':   summary_rio_med,
        'summary_rio_5th':   summary_rio_5th,
        'summary_rio_95th':  summary_rio_95th,
        'tt_med':            tt_med,
        'tt_5th':            tt_5th,
        'tt_95th':           tt_95th,
    }
    print(f'  Done. {n_years} year(s), {len(file_records)} file(s) processed.')

print('\nData extraction complete.')


# =============================================================================
# PLOTTING
# =============================================================================

# Month start positions (0-based day-of-year) and labels for x-axis ticks
_MONTH_DOY   = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
_MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


def make_heatmap(data, years, vessel_class, stat_label, plot_mode,
                 route_label, save_figs, output_dir, total_dist_nm=None,
                 center_doy=258, rio_summary_method='min'):
    """
    Produce a single heatmap figure.

    Parameters
    ----------
    data               : ndarray, shape (n_years, 366)
                         transit_time mode : transit times in days (NaN = impassable)
                         RIO / passability : minimum or summary RIO values
    years              : list of int
    vessel_class       : str
    stat_label         : str   e.g. 'Median', '5th percentile'
    plot_mode          : str   'transit_time', 'RIO', or 'passability'
    route_label        : str   e.g. 'EJCD'
    save_figs          : bool
    output_dir         : str
    total_dist_nm      : float  total route length in nautical miles
    center_doy         : int    0-based DOY to place at the centre of the x-axis
    rio_summary_method : str    'min' or 'mean'; used for RIO colorbar label
    """
    n_years  = len(years)
    n_days   = 365         # show days 0-364; omit leap-year slot at index 365
    plot_arr = data[:, :n_days].copy()

    # Roll columns so that center_doy sits in the middle of the x-axis
    start_doy = (center_doy - n_days // 2) % n_days
    plot_arr  = np.roll(plot_arr, -start_doy, axis=1)

    # Month positions in the rolled frame
    shifted_starts = [(_MONTH_DOY[m] - start_doy) % n_days for m in range(12)]
    month_order    = sorted(range(12), key=lambda m: shifted_starts[m])

    fig_height = max(3.5, n_years * 0.50 + 2.0)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    extent = [-0.5, n_days - 0.5, -0.5, n_years - 0.5]

    # ---------------------------------------------------------------------- #
    # Render according to plot_mode
    # ---------------------------------------------------------------------- #
    if plot_mode == 'transit_time':
        cmap = plt.cm.plasma.copy()
        cmap.set_bad(color='lightgrey')
        masked = np.ma.masked_invalid(plot_arr)
        cfg  = VESSEL_CONFIG[vessel_class]
        vmin = total_dist_nm / cfg['normal_kt'] / 24.0 if total_dist_nm else 0.0
        vmax = 365.0
        im = ax.imshow(masked, aspect='auto', origin='lower',
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=extent, interpolation='none')
        cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
        cbar.set_label('Transit time (days)', fontsize=11)

    elif plot_mode == 'RIO':
        cmap = plt.cm.RdBu.copy()
        cmap.set_bad(color='lightgrey')
        masked  = np.ma.masked_invalid(plot_arr)
        valid   = masked.compressed()
        abs_max = float(np.nanpercentile(np.abs(valid), 98)) if len(valid) > 0 else 30.0
        im = ax.imshow(masked, aspect='auto', origin='lower',
                       cmap=cmap, vmin=-abs_max, vmax=abs_max,
                       extent=extent, interpolation='none')
        cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
        rio_lbl = 'Route-mean RIO' if rio_summary_method == 'mean' else 'Route-minimum RIO (bottleneck)'
        cbar.set_label(rio_lbl, fontsize=11)
        if np.any(~np.isnan(plot_arr)):
            # Contour at RIO = 0 (normal/restricted boundary)
            ax.contour(np.arange(n_days), np.arange(n_years),
                       np.ma.masked_invalid(plot_arr),
                       levels=[0], colors='black', linewidths=0.8, linestyles='--')
            # Contour at RIO = -10 (impassability boundary for PC* vessels)
            ax.contour(np.arange(n_days), np.arange(n_years),
                       np.ma.masked_invalid(plot_arr),
                       levels=[-10], colors='black', linewidths=0.8, linestyles='-')

    elif plot_mode == 'passability':
        thresh = VESSEL_CONFIG[vessel_class]['impass_thresh']
        binary = np.where(
            (~np.isnan(plot_arr)) & (plot_arr >= thresh), 1.0, 0.0
        )
        cmap_bin = mcolors.ListedColormap(['#d73027', '#1a9850'])
        ax.imshow(binary, aspect='auto', origin='lower',
                  cmap=cmap_bin, vmin=0, vmax=1,
                  extent=extent, interpolation='none')
        green_patch = mpatches.Patch(color='#1a9850', label='Passable')
        red_patch   = mpatches.Patch(color='#d73027', label='Impassable')
        ax.legend(handles=[green_patch, red_patch],
                  loc='upper right', fontsize=10, framealpha=0.9)
    else:
        raise ValueError(f"plot_mode must be 'transit_time', 'RIO', or 'passability'; "
                         f"got '{plot_mode}'")

    # ---------------------------------------------------------------------- #
    # Axes formatting
    # ---------------------------------------------------------------------- #
    # x-axis: tick marks, labels, and vertical lines all at the first day of each month
    tick_positions = [shifted_starts[m] for m in month_order]
    tick_labels    = [_MONTH_NAMES[m] for m in month_order]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_xlim(-0.5, n_days - 0.5)

    for m in range(12):
        pos = shifted_starts[m]
        if 0 < pos < n_days:
            ax.axvline(x=pos, color='white', linewidth=0.6, alpha=0.5)

    # y-axis: one tick per year, earliest year at the bottom
    ax.set_yticks(range(n_years))
    ax.set_yticklabels([str(y) for y in years], fontsize=10)
    ax.set_ylabel('Year', fontsize=11)
    ax.set_xlabel('Month', fontsize=11)

    mode_titles = {
        'transit_time': 'Transit Time',
        'RIO':          'RIO',
        'passability':  'Go / No-Go',
    }
    ax.set_title(
        f'{vessel_class}  \u2014  {stat_label}  |  '
        f'{mode_titles[plot_mode]}  |  Route {route_label}',
        fontsize=12, fontweight='bold'
    )

    plt.tight_layout()

    if save_figs:
        os.makedirs(output_dir, exist_ok=True)
        stat_slug = (stat_label.lower()
                     .replace(' ', '_')
                     .replace('th', '')
                     .replace('st', ''))
        fname = os.path.join(
            output_dir,
            f'{vessel_class}_route{route_label}_{plot_mode}_{stat_slug}.png'
        )
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        print(f'  Saved: {fname}')

    plt.close(fig)


# =============================================================================
# GENERATE PLOTS
# =============================================================================

for vessel_class, res in results.items():
    print(f'\nPlotting {vessel_class}  (mode: {plot_mode}) ...')

    if plot_mode == 'transit_time':
        data_trio = [
            (res['tt_med'],           'Median'),
            (res['tt_5th'],           '5th percentile'),
            (res['tt_95th'],          '95th percentile'),
        ]
    elif plot_mode == 'passability':
        # Passability always uses the route-minimum RIO
        data_trio = [
            (res['min_rio_med'],      'Median'),
            (res['min_rio_5th'],      '5th percentile'),
            (res['min_rio_95th'],     '95th percentile'),
        ]
    elif plot_mode == 'RIO':
        # RIO heatmap uses whichever summary method is selected in CONFIG
        data_trio = [
            (res['summary_rio_med'],  'Median'),
            (res['summary_rio_5th'],  '5th percentile'),
            (res['summary_rio_95th'], '95th percentile'),
        ]
    else:
        raise ValueError(f"Unknown plot_mode: '{plot_mode}'")

    for data, stat_label in data_trio:
        make_heatmap(
            data               = data,
            years              = res['years'],
            vessel_class       = vessel_class,
            stat_label         = stat_label,
            plot_mode          = plot_mode,
            route_label        = route_label,
            save_figs          = save_figs,
            output_dir         = output_dir,
            total_dist_nm      = total_dist_nm,
            center_doy         = plot_center_doy,
            rio_summary_method = rio_summary_method,
        )

print('\nDone.')
