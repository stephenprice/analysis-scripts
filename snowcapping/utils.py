"""
Shared utility functions for the E3SM snowcapping analysis scripts.

Provides:
  - Fortran namelist parsing
  - Simulation metadata extraction (latband, dtime, start date)
  - File discovery helpers
  - lnd.log diagnostic parsing
  - ELM history file sorting (handles YYYY-Mon.nc and YYYY-MM-DD-HHHHH.nc)
  - Coordinate normalization helpers
"""

import os
import glob
import gzip
import re

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Fortran namelist parsing
# ---------------------------------------------------------------------------

def parse_namelist(filepath):
    """
    Parse a Fortran namelist file into a nested dict:
        {group_name_lower: {key_lower: raw_value_string}}

    Multi-line values (like hist_fincl lists) store only the first line's
    value; scalar parameters (dtime, start_ymd, etc.) are always correct.
    """
    result = {}
    current_group = None
    try:
        with open(filepath) as fh:
            for raw in fh:
                line = re.sub(r'\s*!.*$', '', raw.strip())   # strip comments
                if not line:
                    continue
                if line.startswith('&'):
                    current_group = line[1:].strip().lower()
                    result[current_group] = {}
                elif line == '/':
                    current_group = None
                elif current_group is not None and '=' in line:
                    key, _, val = line.partition('=')
                    key = key.strip().lower()
                    val = val.strip().rstrip(',').strip()
                    # Keep only the first assignment for each key
                    if key not in result[current_group]:
                        result[current_group][key] = val
    except FileNotFoundError:
        print(f"  Warning: namelist not found: {filepath}")
    return result


def _get_lnd_in(run_dir):
    return parse_namelist(os.path.join(run_dir, 'lnd_in'))


def _get_drv_in(run_dir):
    return parse_namelist(os.path.join(run_dir, 'drv_in'))


# ---------------------------------------------------------------------------
# Simulation metadata
# ---------------------------------------------------------------------------

def get_latband_degrees(run_dir):
    """
    Return the snowcapping latitude band half-width (float) from lnd_in, or
    None if the feature is not active in this run.
    """
    nl = _get_lnd_in(run_dir)
    grp = nl.get('elm_inparm', {})
    active = grp.get('convert_ice_to_river_runoff_latband', '.false.')
    if '.true.' not in active.lower():
        return None
    val = grp.get('convert_ice_to_river_runoff_latband_width_degrees')
    if val is None:
        return None
    try:
        return float(val)
    except ValueError:
        return None


def get_dtime(run_dir):
    """Return the land-model timestep in seconds (int) from lnd_in."""
    nl = _get_lnd_in(run_dir)
    val = nl.get('elm_inparm', {}).get('dtime', '1800')
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return 1800


def get_sim_start_date(run_dir):
    """
    Return the simulation start as a pd.Timestamp, parsed from drv_in
    (seq_timemgr_inparm / start_ymd + start_tod).
    """
    nl = _get_drv_in(run_dir)
    grp = nl.get('seq_timemgr_inparm', {})
    ymd = grp.get('start_ymd', '19000101').strip()
    tod = grp.get('start_tod', '0').strip()
    try:
        ymd = ymd.zfill(8)
        year  = int(ymd[0:4])
        month = int(ymd[4:6])
        day   = int(ymd[6:8])
        sec   = int(float(tod))
        hour, rem = divmod(sec, 3600)
        minute, second = divmod(rem, 60)
        return pd.Timestamp(year=year, month=month, day=day,
                            hour=hour, minute=minute, second=second)
    except Exception:
        print("  Warning: could not parse start_ymd from drv_in; defaulting to 1900-01-01")
        return pd.Timestamp('1900-01-01')


def resolve_latband(config_val, fix_run_dir):
    """
    Return the latitude half-band to use for ice time series metrics.

    Priority:
      1. config_val if not None (user override)
      2. auto-read from fix_run_dir/lnd_in
      3. fall back to 65.0 degrees
    """
    if config_val is not None:
        print(f"  Using user-specified latband: ±{config_val}°")
        return float(config_val)
    val = get_latband_degrees(fix_run_dir)
    if val is not None:
        print(f"  Auto-detected latband from lnd_in: ±{val}°")
        return val
    print("  Warning: could not detect latband from lnd_in; defaulting to ±65°")
    return 65.0


def resolve_map_extent(lat_min, lat_max, lon_min, lon_max, fix_run_dir):
    """
    Return (lat_min, lat_max, lon_min, lon_max) floats for map plots.

    Any None values are filled:
      - lat limits → ±latband from fix run's lnd_in (or ±65°)
      - lon limits → −180 / 180
    User-supplied non-None values are always respected (enables zoom-in).
    """
    if lat_min is None or lat_max is None:
        lb = get_latband_degrees(fix_run_dir)
        lb = lb if lb is not None else 65.0
        lat_min = lat_min if lat_min is not None else -lb
        lat_max = lat_max if lat_max is not None else  lb
    lon_min = lon_min if lon_min is not None else -180.0
    lon_max = lon_max if lon_max is not None else  180.0
    print(f"  Map extent: lat [{lat_min:.1f}, {lat_max:.1f}], "
          f"lon [{lon_min:.1f}, {lon_max:.1f}]")
    return float(lat_min), float(lat_max), float(lon_min), float(lon_max)


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_files(run_dir, pattern, required=True):
    """
    Return a sorted list of files in run_dir matching the glob pattern.
    Prints a warning (does not raise) if no files are found and required=True.
    Works for multi-year runs automatically — all matching files are returned.
    """
    full_pattern = os.path.join(run_dir, pattern)
    files = sorted(glob.glob(full_pattern))
    if not files and required:
        print(f"  Warning: no files matching '{full_pattern}'")
    return files


# ---------------------------------------------------------------------------
# lnd.log diagnostic parsing
# ---------------------------------------------------------------------------

# Compiled regexes for per-step and cumulative diagnostic lines.
# Format from SnowHydrologyMod.F90 SnowCappingDiagLog():
#   write(iulog,'(a,i10,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,1pe12.4,a,i10)')
#       'SNOWCAP_LATBAND_DIAG step=', nstep,
#       ' step_ice_mass_kg=', ..., ' step_latent_energy_j=', ...,
#       ' step_mean_cooling_k=', ..., ' step_max_cooling_k=', ...,
#       ' step_cols=', g_step_cols
#
# The 1pe12.4 Fortran format produces values like " 1.2345e+06" (positive,
# leading space) or "-1.2345e-02" (negative, no leading space).

_FP = r'\s*([-]?[\d.]+[eE][+-]?\d+)'   # one Fortran 1pe12.4 value

_STEP_RE = re.compile(
    r'SNOWCAP_LATBAND_DIAG step=\s*(\d+)'
    r'\s+step_ice_mass_kg='    + _FP +
    r'\s+step_latent_energy_j=' + _FP +
    r'\s+step_mean_cooling_k='  + _FP +
    r'\s+step_max_cooling_k='   + _FP +
    r'\s+step_cols=\s*(\d+)'
)

_CUM_RE = re.compile(
    r'SNOWCAP_LATBAND_DIAG_CUM step=\s*(\d+)'
    r'\s+cum_ice_mass_kg='     + _FP +
    r'\s+cum_latent_energy_j=' + _FP +
    r'\s+cum_mean_cooling_k='  + _FP +
    r'\s+cum_max_cooling_k='   + _FP +
    r'\s+cum_cols=\s*(\d+)'
)


def parse_lnd_logs(run_dir, start_date=None, dtime_seconds=1800):
    """
    Parse all lnd.log.*.gz files in run_dir for SNOWCAP_LATBAND_DIAG entries.

    Parameters
    ----------
    run_dir        : path to the simulation run directory
    start_date     : pd.Timestamp for simulation start (read from drv_in if None)
    dtime_seconds  : land-model timestep in seconds (read from lnd_in if not supplied)

    Returns
    -------
    step_df, cum_df : pd.DataFrames indexed by model datetime.
        Both are empty DataFrames if no diagnostic lines are found.

    Notes
    -----
    - Each lnd.log.*.gz file covers one job segment.  Files are sorted by
      filename (job IDs are monotonically increasing) so step order is correct.
    - The model datetime is computed as:
          start_date + timedelta(seconds = nstep * dtime_seconds)
      where nstep is the cumulative step count from the start of the simulation.
    """
    if start_date is None:
        start_date = get_sim_start_date(run_dir)

    log_files = sorted(glob.glob(os.path.join(run_dir, 'lnd.log.*.gz')))
    if not log_files:
        print(f"  No lnd.log.*.gz files found in {run_dir}")
        return pd.DataFrame(), pd.DataFrame()

    step_rows, cum_rows = [], []

    for logfile in log_files:
        try:
            with gzip.open(logfile, 'rt', encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    # Check per-step line first (more common)
                    m = _STEP_RE.search(line)
                    if m:
                        nstep = int(m.group(1))
                        dt = start_date + pd.Timedelta(seconds=nstep * dtime_seconds)
                        step_rows.append({
                            'datetime':           dt,
                            'step':               nstep,
                            'step_ice_mass_kg':   float(m.group(2)),
                            'step_latent_energy_j': float(m.group(3)),
                            'step_mean_cooling_k':  float(m.group(4)),
                            'step_max_cooling_k':   float(m.group(5)),
                            'step_cols':            int(m.group(6)),
                        })
                        continue
                    m = _CUM_RE.search(line)
                    if m:
                        nstep = int(m.group(1))
                        dt = start_date + pd.Timedelta(seconds=nstep * dtime_seconds)
                        cum_rows.append({
                            'datetime':          dt,
                            'step':              nstep,
                            'cum_ice_mass_kg':   float(m.group(2)),
                            'cum_latent_energy_j': float(m.group(3)),
                            'cum_mean_cooling_k':  float(m.group(4)),
                            'cum_max_cooling_k':   float(m.group(5)),
                            'cum_cols':            int(m.group(6)),
                        })
        except Exception as exc:
            print(f"  Warning: error reading {os.path.basename(logfile)}: {exc}")

    step_df = (pd.DataFrame(step_rows).set_index('datetime')
               if step_rows else pd.DataFrame())
    cum_df  = (pd.DataFrame(cum_rows).set_index('datetime')
               if cum_rows  else pd.DataFrame())

    if step_df.empty:
        print(f"  No SNOWCAP_LATBAND_DIAG entries found in {run_dir}")
    else:
        print(f"  Parsed {len(step_df)} SNOWCAP_LATBAND_DIAG entries "
              f"({len(log_files)} log files) from {run_dir}")

    return step_df, cum_df


# ---------------------------------------------------------------------------
# ELM history file helpers
# ---------------------------------------------------------------------------

# Standard ELM monthly history filename:
#   *.elm.hN.YYYY-MM-DD-HHHHH.nc
# The timestamp marks the *start* of the next month, so YYYY-02-01-00000
# contains the January monthly average.  Month = MM - 1 with year wrap.
_ELM_TS_RE = re.compile(r'\.elm\.h\d+\.(\d{4})-(\d{2})-\d{2}-\d{5}\.nc$')


def month_from_elm_filename(fname):
    """
    Extract (year, month) integers from a standard ELM monthly history file.

    Accepted pattern:
      *.elm.hN.YYYY-MM-DD-HHHHH.nc  (end-of-averaging-period timestamp)
        e.g. elm.h0.1900-02-01-00000.nc  → January 1900 average

    Returns (year, month) or (None, None) if the pattern is not recognised.
    Non-standard files such as YYYY-Mon.nc or YYYY-July.nc are intentionally
    not matched here so they are silently skipped by get_elm_files_by_month.
    """
    m = _ELM_TS_RE.search(os.path.basename(fname))
    if not m:
        return None, None
    year, mm = int(m.group(1)), int(m.group(2))
    # Shift back one month: the file written at YYYY-MM-01 holds month MM-1
    mm -= 1
    if mm == 0:
        mm, year = 12, year - 1
    return year, mm


def get_elm_files_by_month(run_dir, stream):
    """
    Discover standard ELM monthly history files for a given output stream
    (e.g. 'h0', 'h2'), excluding restart-history files (elm.rh0.*, …) and
    any non-standard derived files (e.g. YYYY-Mon.nc, YYYY-July.nc).

    Returns a dict:  (year, month) -> filepath
    Works for multi-year runs (returns all years present in run_dir).
    """
    all_files = find_files(run_dir, f'*.elm.{stream}.*.nc')
    # Exclude restart-history files: .elm.rh0.*, .elm.rh2.*, etc.
    all_files = [f for f in all_files
                 if f'.elm.r{stream}.' not in os.path.basename(f)]

    by_month = {}
    for fpath in all_files:
        year, month = month_from_elm_filename(fpath)
        if year is None:
            # Non-standard file (e.g. YYYY-Mon.nc) — skip silently
            continue
        key = (year, month)
        if key not in by_month:
            by_month[key] = fpath

    if not by_month:
        print(f"  Warning: no elm.{stream} history files found in {run_dir}")
    else:
        print(f"  Found {len(by_month)} elm.{stream} file(s) in {run_dir}")

    return by_month


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def normalize_lon(lon):
    """Convert longitude array to the range [−180, 180]."""
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def rad_to_deg_if_needed(arr):
    """
    Convert an angular array from radians to degrees if the values suggest
    radian units.  Two conventions are handled:

    1. Symmetric range [-π, π] (or sub-range thereof, e.g. latitude): detected
       when max(|arr|) ≤ π + ε.
    2. Unsigned [0, 2π] range — the MPAS convention for lonCell: detected when
       min(arr) ≥ -ε and max(arr) ≤ 2π + ε.
    """
    arr = np.asarray(arr, dtype=float)
    amax = np.nanmax(np.abs(arr))
    if amax <= np.pi + 0.02:
        return np.degrees(arr)
    # MPAS lonCell: stored as unsigned radians in [0, 2π]
    if np.nanmin(arr) >= -0.02 and amax <= 2 * np.pi + 0.02:
        return np.degrees(arr)
    return arr


def parse_mpas_xtime(ds):
    """
    Parse an MPAS xtime character variable from an open xarray Dataset into
    a list of pd.Timestamps.

    Tries several candidate variable names used by different MPAS analysis
    members (timeSeriesStatsDaily, timeSeriesStatsMonthly, etc.).
    Returns a list of pd.Timestamp, or None on failure.
    """
    _XTIME_CANDS = [
        'xtime_startDaily', 'xtime_endDaily',
        'xtime_startMonthly', 'xtime_endMonthly',
        'xtime_start', 'xtime_end', 'xtime',
    ]
    raw = None
    for vname in _XTIME_CANDS:
        if vname in ds:
            raw = ds[vname].values
            break
    if raw is None:
        return None
    try:
        # shape (Time,) of bytes/str  or  (Time, StrLen) of individual chars
        if raw.ndim == 2:
            strs = [
                ''.join(c.decode('utf-8') if isinstance(c, bytes) else str(c)
                        for c in row).strip()
                for row in raw
            ]
        else:
            strs = [
                s.decode('utf-8').strip() if isinstance(s, bytes) else str(s).strip()
                for s in raw
            ]
        # Format: "YYYY-MM-DD_HH:MM:SS" → replace underscore, take first 19 chars
        strs = [s[:19].replace('_', ' ') for s in strs]
        return list(pd.to_datetime(strs))
    except Exception as exc:
        print(f"  Warning: could not parse MPAS xtime ({exc})")
        return None


# ---------------------------------------------------------------------------
# MPAS mesh coordinate loading
# ---------------------------------------------------------------------------

def load_mpassi_mesh(run_dir):
    """
    Load latCell, lonCell, and areaCell from an mpassi restart file.

    These mesh variables are not written to the analysis-member history files
    (timeSeriesStatsDaily, timeSeriesStatsMonthly) to save disk space, but
    are always present in the restart files (*.mpassi.rst.*.nc).

    Returns (lat_deg, lon_deg, area_m2) as 1-D numpy arrays, or raises
    RuntimeError if no restart file can be found or opened.
    """
    rst_files = find_files(run_dir, '*.mpassi.rst.*.nc', required=False)
    # Prefer the earliest restart (smallest step count = closest to initial mesh)
    rst_files = [f for f in rst_files if '.rst.am.' not in f]  # skip analysis-member restarts
    if not rst_files:
        raise RuntimeError(
            f"No mpassi restart file found in {run_dir}.\n"
            "Cannot determine MPAS sea-ice mesh coordinates (latCell, areaCell)."
        )

    rst_file = rst_files[0]
    print(f"    Reading MPASSI mesh from: {os.path.basename(rst_file)}")
    try:
        import xarray as xr
        ds = xr.open_dataset(rst_file, decode_times=False, mask_and_scale=False)
    except Exception as exc:
        raise RuntimeError(f"Cannot open mpassi restart file: {exc}")

    # Latitude
    lat = None
    for v in ['latCell', 'lat']:
        if v in ds:
            lat = rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lat is None:
        raise RuntimeError(
            f"No latCell variable in {os.path.basename(rst_file)}. "
            f"Available: {list(ds.data_vars)[:20]}"
        )

    # Longitude
    lon = None
    for v in ['lonCell', 'lon']:
        if v in ds:
            lon = rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lon is None:
        raise RuntimeError(f"No lonCell variable in {os.path.basename(rst_file)}.")
    lon = normalize_lon(lon)

    # Cell area (m²)
    area = None
    for v in ['areaCell', 'area']:
        if v in ds:
            area = ds[v].values.ravel().astype(float)
            break
    if area is None:
        print(f"  Warning: no areaCell in restart file; using uniform weights.")
        area = np.ones(len(lat), dtype=float)

    ds.close()
    print(f"    Mesh: {len(lat)} cells")
    return lat, lon, area


def load_mpaso_mesh(run_dir):
    """
    Load latCell, lonCell, and areaCell from an mpaso restart file.

    Returns (lat_deg, lon_deg, area_m2) as 1-D numpy arrays.
    """
    rst_files = find_files(run_dir, '*.mpaso.rst.*.nc', required=False)
    rst_files = [f for f in rst_files if '.rst.am.' not in f]
    if not rst_files:
        raise RuntimeError(
            f"No mpaso restart file found in {run_dir}.\n"
            "Cannot determine MPAS ocean mesh coordinates."
        )

    rst_file = rst_files[0]
    print(f"    Reading MPASO mesh from: {os.path.basename(rst_file)}")
    try:
        import xarray as xr
        ds = xr.open_dataset(rst_file, decode_times=False, mask_and_scale=False)
    except Exception as exc:
        raise RuntimeError(f"Cannot open mpaso restart file: {exc}")

    lat = None
    for v in ['latCell', 'lat']:
        if v in ds:
            lat = rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lat is None:
        raise RuntimeError(f"No latCell in {os.path.basename(rst_file)}.")

    lon = None
    for v in ['lonCell', 'lon']:
        if v in ds:
            lon = rad_to_deg_if_needed(ds[v].values.ravel())
            break
    if lon is None:
        raise RuntimeError(f"No lonCell in {os.path.basename(rst_file)}.")
    lon = normalize_lon(lon)

    area = None
    for v in ['areaCell', 'area']:
        if v in ds:
            area = ds[v].values.ravel().astype(float)
            break
    if area is None:
        print(f"  Warning: no areaCell in mpaso restart file; using uniform weights.")
        area = np.ones(len(lat), dtype=float)

    ds.close()
    print(f"    Mesh: {len(lat)} cells")
    return lat, lon, area
