"""
E3SM Snowcapping Fix — Analysis Configuration
==============================================

Edit the settings in this file before running any of the plot scripts.

Usage
-----
1.  Load the E3SM unified Python environment:
        source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_login.sh

2.  Edit FIX_RUN_DIR / CTRL_RUN_DIR (and optionally the other settings below).

3.  Run any of the three plotting scripts independently:
        python plot_timeseries_ice.py      # daily sea-ice metrics time series
        python plot_timeseries_lndlog.py   # SNOWCAP_LATBAND_DIAG from lnd.log
        python plot_maps.py                # monthly map-view for all 5 variables

Notes
-----
- All three scripts import this file, so changes here apply everywhere.
- The scripts work for multi-year runs automatically (they glob all matching files).
- Set any of the None values below to override auto-detection from the namelists.
"""

# ---------------------------------------------------------------------------
# Simulation run directories
# Point each to the 'run/' subdirectory of the respective E3SM case.
# ---------------------------------------------------------------------------

FIX_RUN_DIR = (
    "/pscratch/sd/s/sprice/e3sm_scratch/pm-cpu/"
    "20260427.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40.pm-cpu.snowcappingFix/run/"
)

CTRL_RUN_DIR = (
    "/pscratch/sd/s/sprice/e3sm_scratch/pm-cpu/"
    "20260427.BGWCYCL2010.ne30pg2_r05_IcoswISC30E3r5_gis4to40.pm-cpu.snowcappingFixBaseline/run/"
)

# ---------------------------------------------------------------------------
# Labels used in plot legends and titles
# ---------------------------------------------------------------------------

FIX_LABEL  = "SnowcappingFix (±65°)"
CTRL_LABEL = "Control (no fix)"

# ---------------------------------------------------------------------------
# Latitude band for sea-ice time series
#
# Set to None  → auto-read convert_ice_to_river_runoff_latband_width_degrees
#                from FIX_RUN_DIR/lnd_in (recommended).
# Set to float → override, e.g. LATBAND_DEGREES = 45.0 to restrict analysis
#                to ±45° regardless of what the namelist says.
#
# The same band is applied to BOTH simulations so the comparison is fair.
# ---------------------------------------------------------------------------

LATBAND_DEGREES = None  # e.g. 65.0 or None
#LATBAND_DEGREES = 35.0  # e.g. 65.0 or None

# ---------------------------------------------------------------------------
# Map plot spatial extent (degrees)
#
# Set any value to None → auto-derive from the fix run's lnd_in latband (lat)
#                          or default to ±180° (lon).
# Override any or all to zoom into a sub-region, e.g.:
#   MAP_LAT_MIN = 40.0
#   MAP_LAT_MAX = 45.0
#   MAP_LON_MIN = -130.0
#   MAP_LON_MAX = -60.0
# ---------------------------------------------------------------------------

MAP_LAT_MIN = None   # southern limit  (degrees_north)
MAP_LAT_MAX = None   # northern limit  (degrees_north)
MAP_LON_MIN = None   # western  limit  (degrees_east, range −180 to 180)
MAP_LON_MAX = None   # eastern  limit  (degrees_east)

# test focus area: zoom to high mountain Asia
#MAP_LAT_MIN = 0.0    # southern limit  (degrees_north)
#MAP_LAT_MAX = 50.0   # northern limit  (degrees_north)
#MAP_LON_MIN = 45.0   # western  limit  (degrees_east, range −180 to 180)
#MAP_LON_MAX = 105.0  # eastern  limit  (degrees_east)

# test focus area: zoom to tip of S. Greenland 
#MAP_LAT_MIN = 50.0    # southern limit  (degrees_north)
#MAP_LAT_MAX = 70.0    # northern limit  (degrees_north)
#MAP_LON_MIN = -60.0   # western  limit  (degrees_east, range −180 to 180)
#MAP_LON_MAX = -30.0   # eastern  limit  (degrees_east)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUTPUT_DIR = "./plots"   # directory for all saved figures (created if absent)

# ---------------------------------------------------------------------------
# Plot styling
# ---------------------------------------------------------------------------

COLORS     = {"fix": "steelblue",  "ctrl": "firebrick"}
LINESTYLES = {"fix": "-",          "ctrl": "--"}
LINEWIDTHS = {"fix": 1.5,          "ctrl": 1.5}

DPI = 150   # figure resolution for saved PNGs

# ---------------------------------------------------------------------------
# Per-variable colorbar overrides for map plots
#
# Set entries here to manually clamp the colour range for specific variables.
# Keys must match the 'name' field in MAP_VARIABLES (plot_maps.py).
# Values are (vmin, vmax) tuples.  Set to None to use auto-scaling.
#
# Example:
#   MAP_VRANGE_OVERRIDES = {
#       'iceRunoffFlux': (0.0, 3e-4),
#       'riverRunoffFlux': (0.0, 1e-3),
#   }
# ---------------------------------------------------------------------------

#MAP_VRANGE_OVERRIDES = {}   # empty → use auto-scaling for all variables
MAP_VRANGE_OVERRIDES = {
    'iceRunoffFlux': (0.0, 4e-24),
}

# ---------------------------------------------------------------------------
# Select which map variables to plot
#
# Set to None → plot all variables (default).
# Set to a list of variable names to plot only those, e.g.:
#   MAP_VARIABLES_TO_PLOT = ['iceAreaCell', 'QSNWCPICE']
#
# Valid names: iceAreaCell, iceRunoffFlux, riverRunoffFlux, QSNWCPICE, SNO_T_davg
# ---------------------------------------------------------------------------

MAP_VARIABLES_TO_PLOT = None   # None → plot all
#MAP_VARIABLES_TO_PLOT = [ 'iceAreaCell' ]

# ---------------------------------------------------------------------------
# Sea-ice concentration masking threshold for map plots
#
# Values of iceAreaCell at or below this threshold are masked out (rendered as
# transparent background) so the colorbar only colours cells with meaningful
# sea-ice concentration.  This makes it much easier to distinguish very low
# concentrations from truly ice-free ocean.
#
# Set to 0.0 (default) to mask only identically-zero cells.
# Set to a small positive value (e.g. 1e-6) to also exclude near-zero noise.
# Set to None to disable masking entirely (original behaviour).
# ---------------------------------------------------------------------------

ICE_MASK_THRESHOLD = 0.0
