"""
Shared helpers for climate-year (Mar 15 -> Mar 14) trend analysis.

Used by:
    RIO_trend_climate_year.ipynb
    seaice_trend_climate_year.ipynb

Both notebooks load a "climate-year" NetCDF dataset (produced by the export
cell added to RIO_daily_heatmaps.ipynb / seaice_daily_heatmaps.ipynb) in which
DOY 0 = March 15 and DOY 364 = March 14 of the following calendar year. Because
every season is fully contained within a single 365-day row, none of the
Dec/Jan year-boundary wraparound handling that was needed for calendar-year
data is required here -- `find_crossing_doy` below has no wrap parameter at
all, by construction.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import linregress as _linregress

# ---------------------------------------------------------------------------
# Shifted-DOY month ticks (DOY 0 = March 15)
# ---------------------------------------------------------------------------
# Calendar month-start DOYs (0-based, non-leap year): Jan=0, Feb=31, Mar=59,
# Apr=90, May=120, Jun=151, Jul=181, Aug=212, Sep=243, Oct=273, Nov=304, Dec=334.
# Shifted = (calendar_doy - 73) % 365. Listed in chronological order starting
# from the first full month after March 15 (April).
DOY_SHIFT = 73  # 0-based calendar DOY of March 15 (non-leap year)

SHIFTED_MONTH_DOY = [17, 47, 78, 108, 139, 170, 200, 231, 261, 292, 323, 351]
SHIFTED_MONTH_NAMES = [
    'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
    'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar',
]


def calendar_doy_to_shifted(calendar_doy, shift=DOY_SHIFT):
    """
    Translate a familiar 0-based calendar day-of-year into shifted (climate-year)
    DOY space, where 0 = March 15.

    Parameters
    ----------
    calendar_doy : int or array-like   0-based calendar DOY (0 = Jan 1)
    shift        : int                 defaults to DOY_SHIFT (73 = March 15)

    Returns
    -------
    int or ndarray : shifted DOY, always in [0, 364]
    """
    return (np.asarray(calendar_doy) - shift) % 365


# ---------------------------------------------------------------------------
# Smoothing, crossing-detection, and trend-fitting helpers
# ---------------------------------------------------------------------------
def smooth_by_year(data, window):
    """
    Apply a centred rolling average along the DOY axis for each season.

    NaN values (missing days) are excluded from each window mean. Edge days use
    a one-sided window rather than reflection padding so that NaN values at the
    end of the array are not mirrored inward.

    Parameters
    ----------
    data   : ndarray, shape (n_seasons, n_doys)
    window : int   number of days (odd recommended)

    Returns
    -------
    smoothed : ndarray, same shape, NaN where the window is entirely NaN
    """
    if window <= 1:
        return data.copy()
    half = window // 2
    n_seasons, n_cols = data.shape
    smoothed = np.full_like(data, np.nan, dtype=float)
    for yi in range(n_seasons):
        row = data[yi].astype(float)
        for d in range(n_cols):
            d0 = max(0, d - half)
            d1 = min(n_cols, d + half + 1)
            vals = row[d0:d1]
            ok = vals[~np.isnan(vals)]
            smoothed[yi, d] = np.mean(ok) if len(ok) > 0 else np.nan
    return smoothed


def find_crossing_doy(smoothed_row, threshold, window_start, window_end, find):
    """
    Find the DOY at which a smoothed time-series crosses a threshold.

    Unlike the calendar-year version of this function, there is no year-boundary
    wraparound to handle here: every season is fully self-contained in one row
    of shifted-DOY space (0 = March 15 .. 364 = March 14 of the following year),
    so `window_end` never needs to exceed `len(smoothed_row) - 1`.

    Parameters
    ----------
    smoothed_row  : ndarray (n_days,)  smoothed daily values for one season
    threshold     : float
    window_start  : int   first DOY to search (0-based, inclusive)
    window_end    : int   last  DOY to search (0-based, inclusive)
    find          : str
        'last_above'  - last  DOY where smoothed_row >= threshold
        'first_above' - first DOY where smoothed_row >= threshold

    Returns
    -------
    int or float(NaN) : crossing DOY, or NaN if the threshold is not crossed
    """
    doys = np.arange(window_start, min(window_end + 1, len(smoothed_row)))
    vals = smoothed_row[doys]
    above = doys[vals >= threshold]
    if len(above) == 0:
        return np.nan
    return int(above[-1]) if find == 'last_above' else int(above[0])


def compute_trend(doys, years):
    """
    Fit a linear trend (DOY ~ season year) and return diagnostics.

    Years with NaN crossing dates are excluded; a warning is printed listing them.

    Parameters
    ----------
    doys  : array-like (n_seasons,)  crossing DOY per season; NaN = no crossing
    years : array-like (n_seasons,)  corresponding season-start calendar years

    Returns
    -------
    slope, intercept : float   trend (days / year) and intercept
    r2               : float   coefficient of determination (R^2)
    rmse             : float   root-mean-square error of residuals (days)
    valid_years      : ndarray years included in the fit
    valid_doys       : ndarray crossing DOYs included in the fit
    """
    doys = np.asarray(doys, dtype=float)
    years = np.asarray(years, dtype=float)
    valid = ~np.isnan(doys)
    excluded = years[~valid].astype(int).tolist()
    if excluded:
        print(f'    WARNING: threshold not crossed in {len(excluded)} season(s): '
              f'{excluded}  -- excluded from trend fit.')
    valid_years = years[valid]
    valid_doys = doys[valid]
    if len(valid_years) < 2:
        print('    WARNING: fewer than 2 valid seasons -- cannot compute trend.')
        return np.nan, np.nan, np.nan, np.nan, valid_years, valid_doys
    res = _linregress(valid_years, valid_doys)
    r2 = res.rvalue ** 2
    predicted = res.slope * valid_years + res.intercept
    rmse = float(np.sqrt(np.mean((valid_doys - predicted) ** 2)))
    return res.slope, res.intercept, r2, rmse, valid_years, valid_doys


# ---------------------------------------------------------------------------
# Generic heatmap + trend-overlay plot
# ---------------------------------------------------------------------------
def plot_climate_year_heatmap_with_trend(
    data, season_years, stat_label, field_label,
    cmap_name='RdBu', vmin=-30, vmax=30, cbar_label='',
    center_doy=185, contour_levels=None,
    trend_overlay=None, route_label='', save_figs=False,
    output_dir='.', field_slug='field',
):
    """
    Produce a single climate-year heatmap figure (one row per season, shifted-DOY
    x-axis) with optional threshold-crossing scatter points and trend lines.

    Because each row represents exactly one continuous ice season, a fall or
    spring crossing lands on exactly one row at exactly one column -- there is
    no possibility of the double/missing-marker ambiguity that calendar-year
    rows produce for wrapped crossings.

    Parameters
    ----------
    data           : ndarray, shape (n_seasons, 365)  already in shifted-DOY space
    season_years   : list of int   season-start calendar year per row
    stat_label     : str   e.g. 'Median', '5th percentile'
    field_label    : str   e.g. 'RIO', 'Ice Concentration', 'Ice Thickness'
    cmap_name      : str   matplotlib colormap name
    vmin, vmax     : float colorbar limits (vmax=None to auto-derive from data)
    cbar_label     : str   colorbar axis label
    center_doy     : int   shifted DOY (0=Mar15) to place at the centre of the x-axis
    contour_levels : list of float or None   optional contour lines (e.g. [0] for RIO)
    trend_overlay  : dict or None
                     Keys 'spring' and/or 'fall', each a dict with:
                       'yi'        - ndarray of season-row indices with valid crossings
                       'rolled_x'  - corresponding rolled x-axis positions
                       'slope'     - trend slope (days / year)
                       'intercept' - trend intercept
    route_label    : str   used in title/filename
    save_figs      : bool
    output_dir     : str
    field_slug     : str   short identifier for filenames (e.g. 'rio', 'conc', 'thk')

    Returns
    -------
    fig, ax : the created matplotlib Figure and Axes
    """
    import os as _os
    import matplotlib.pyplot as plt
    import numpy as np

    n_seasons = len(season_years)
    n_days = 365
    plot_arr = data[:, :n_days].copy()

    start_doy = (center_doy - n_days // 2) % n_days
    plot_arr = np.roll(plot_arr, -start_doy, axis=1)

    shifted_starts = [(d - start_doy) % n_days for d in SHIFTED_MONTH_DOY]
    month_order = sorted(range(12), key=lambda m: shifted_starts[m])

    if vmax is None:
        vmax = float(np.nanmax(plot_arr)) if np.any(np.isfinite(plot_arr)) else 1.0

    fig_height = max(3.5, n_seasons * 0.50 + 2.0)
    fig, ax = plt.subplots(figsize=(16, fig_height))

    extent = [-0.5, n_days - 0.5, -0.5, n_seasons - 0.5]
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color='lightgrey')
    masked = np.ma.masked_invalid(plot_arr)

    im = ax.imshow(masked, aspect='auto', origin='lower',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    extent=extent, interpolation='none')
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.03)
    cbar.set_label(cbar_label, fontsize=11)

    if contour_levels and np.any(~np.isnan(plot_arr)):
        # matplotlib requires levels passed in a single contour() call to be
        # strictly increasing, and callers may reasonably pass levels in any
        # order (e.g. [0, -10]) with different intended linestyles per level
        # (dashed for the primary threshold, solid for a secondary one), so
        # each level is drawn with its own contour() call rather than one
        # combined call.
        for _i, _level in enumerate(contour_levels):
            ax.contour(np.arange(n_days), np.arange(n_seasons), masked,
                       levels=[_level], colors='black',
                       linewidths=0.8, linestyles=('--' if _i == 0 else '-'))

    tick_positions = [shifted_starts[m] for m in month_order]
    tick_labels = [SHIFTED_MONTH_NAMES[m] for m in month_order]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_xlim(-0.5, n_days - 0.5)

    for m in range(12):
        pos = shifted_starts[m]
        if 0 < pos < n_days:
            ax.axvline(x=pos, color='white', linewidth=0.6, alpha=0.5)

    ax.set_yticks(range(n_seasons))
    ax.set_yticklabels([f'{y}-{y + 1}' for y in season_years], fontsize=10)
    ax.set_ylabel('Season (Mar\u2013Mar)', fontsize=11)
    ax.set_xlabel('Month', fontsize=11)
    ax.set_title(
        f'{field_label}  \u2014  {stat_label}  |  Route {route_label}\n'
        f'Climate year: March 15 \u2013 March 14',
        fontsize=12, fontweight='bold',
    )

    if trend_overlay is not None:
        _overlay_styles = {
            'spring': {'color': 'limegreen', 'marker': 'o', 'label': 'Spring crossing'},
            'fall':   {'color': 'orange',    'marker': 's', 'label': 'Fall crossing'},
        }
        legend_handles = []
        year_arr = np.array(season_years, dtype=float)
        for season_key, style in _overlay_styles.items():
            ov = trend_overlay.get(season_key)
            if ov is None or len(ov.get('yi', [])) == 0:
                continue
            sc = ax.scatter(ov['rolled_x'], ov['yi'],
                            color=style['color'], marker=style['marker'],
                            s=40, zorder=5, edgecolors='black', linewidths=0.5,
                            label=style['label'])
            legend_handles.append(sc)
            if not np.isnan(ov['slope']):
                pred_doy = ov['slope'] * year_arr + ov['intercept']
                pred_rolled = np.array([(d - start_doy) % n_days for d in pred_doy])
                trend_yi = np.arange(n_seasons, dtype=float)
                ax.plot(pred_rolled, trend_yi,
                       color=style['color'], linewidth=2, linestyle='--', zorder=4)
        if legend_handles:
            ax.legend(handles=legend_handles, fontsize=8,
                     loc='lower right', framealpha=0.85)

    plt.tight_layout()

    if save_figs:
        _os.makedirs(output_dir, exist_ok=True)
        stat_slug = (stat_label.lower().replace(' ', '_')
                     .replace('th', '').replace('st', ''))
        overlay_suffix = '_trend' if trend_overlay is not None else ''
        fname = _os.path.join(
            output_dir,
            f'{field_slug}_climate_year_route{route_label}_{stat_slug}{overlay_suffix}.png',
        )
        fig.savefig(fname, dpi=150, bbox_inches='tight')
        print(f'  Saved: {fname}')

    return fig, ax
