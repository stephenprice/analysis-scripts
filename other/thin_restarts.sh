#!/usr/bin/env bash
# thin_restarts.sh
# Thins E3SM monthly restart files to one retained month per year
# within a user-specified year range.

set -euo pipefail

# ----------------------------------------------------------------
# Help / usage
# ----------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $(basename "$0") [--help]

Interactively thins E3SM monthly restart files to one retained month
per year within a user-specified year range.

All inputs are prompted interactively:
  - Restart file directory
  - Whether to search subdirectories recursively
  - Start and end year (YYYY)
  - Month to retain (1–12; all other months in range are candidates for moving)
  - Dry-run mode (lists files to move/keep without touching anything)
  
Recognised E3SM components (detected from filename tokens):
  eam      — atmosphere (EAM)
  elm      — land model (ELM)
  mosart   — river runoff (MOSART)
  mpaso    — ocean (MPAS-Ocean)
  mpassi   — sea ice (MPAS-Seaice)
  mali     — land ice (MALI / MPAS-Albany-LandIce)
  cpl      — coupler
  datm dlnd drof — data components
  
File types handled (restart files only — history/initial files are ignored):
  CIME-style : case.comp.r.YYYY-MM-DD-SSSSS.nc
                case.comp.rs.YYYY-MM-DD-SSSSS.nc
                case.comp.rh0.YYYY-MM-DD-SSSSS.nc  (and rh1, rh2, …)
  MPAS-style : case.comp.rst.YYYY-MM-DD_HH.MM.SS.nc
  
Explicitly skipped (not restart files):
  case.comp.h0.…nc, case.comp.hi.…nc, case.comp.hist.…nc  (history)
  case.comp.i.…nc  (initial conditions)
  
Dry-run output:
  Lists files that WOULD be moved to ./restarts/ (with sizes) and
  files that WOULD be retained, plus a total storage estimate.
  
Examples: 
  $(basename "$0")           # run interactively
  $(basename "$0") --help    # show this message
EOF
  exit 0
} 
  
[[ "${1-}" == "--help" || "${1-}" == "-h" ]] && usage
    
# ----------------------------------------------------------------
# Helper: format bytes into human-readable string
# ----------------------------------------------------------------
format_bytes() {
  local bytes=$1
  if   (( bytes >= 1099511627776 )); then
    awk "BEGIN {printf \"%.2f TB\", $bytes/1099511627776}"
  elif (( bytes >= 1073741824 )); then
    awk "BEGIN {printf \"%.2f GB\", $bytes/1073741824}"
  elif (( bytes >= 1048576 )); then
    awk "BEGIN {printf \"%.2f MB\", $bytes/1048576}"
  else
    echo "${bytes} B"
  fi
}

# ----------------------------------------------------------------
# Helper: extract component name from an E3SM restart filename.
# Matches known component tokens that appear between dots (or at
# the start of the filename for MPAS files without a case prefix).
# ----------------------------------------------------------------
get_component() {
  local fname="$1"
  local known=(eam elm mosart mpaso mpassi mali cpl datm dlnd drof)
  for comp in "${known[@]}"; do
    if [[ "$fname" =~ \."$comp"\. || "$fname" =~ \."$comp"_ || "$fname" == "$comp".* ]]; then
      echo "$comp"
      return
    fi
  done
  echo "unknown"
}

# ----------------------------------------------------------------
# Helper: return 0 (true) if the file is a restart file.
# The token immediately after .<comp>. must be r, rs, rst, or rh<N>.
# This excludes history files (.h0., .hi., .hist., .hist.am. …) and
# initial-conditions files (.i.).
# ----------------------------------------------------------------
is_restart_file() {
  local fname="$1"
  local comp="$2"
  [[ "$fname" =~ \."$comp"\.(r|rs|rst|rh[0-9]+)\. ]]
}

# ================================================================
# Prompt for inputs
# ================================================================
echo "=================================================="
echo "  E3SM Restart File Thinning Script"
echo "=================================================="
echo

read -rp "Restart file directory          : " RST_DIR
RST_DIR="${RST_DIR%/}"
if [[ ! -d "$RST_DIR" ]]; then
  echo "ERROR: Directory not found: $RST_DIR" >&2; exit 1
fi

read -rp "Search subdirectories? [y/N]    : " RECURSE_INPUT
RECURSE=false
[[ "${RECURSE_INPUT,,}" =~ ^y ]] && RECURSE=true

echo
read -rp "Start year (YYYY)               : " START_YEAR
read -rp "End year   (YYYY)               : " END_YEAR
read -rp "Month to RETAIN (1-12)          : " RETAIN_MONTH_RAW

echo
read -rp "Dry run only? [y/N]             : " DRY_RUN_INPUT
DRY_RUN=false
[[ "${DRY_RUN_INPUT,,}" =~ ^y ]] && DRY_RUN=true

# ----------------------------------------------------------------
# Validate inputs
# ----------------------------------------------------------------
if ! [[ "$START_YEAR" =~ ^[0-9]{4}$ && "$END_YEAR" =~ ^[0-9]{4}$ ]]; then
  echo "ERROR: Years must be 4-digit integers." >&2; exit 1
fi
if (( START_YEAR > END_YEAR )); then
  echo "ERROR: Start year must be <= end year." >&2; exit 1
fi
if ! [[ "$RETAIN_MONTH_RAW" =~ ^([1-9]|1[0-2])$ ]]; then
  echo "ERROR: Month must be between 1 and 12." >&2; exit 1
fi

RETAIN_MONTH=$(printf '%02d' "$RETAIN_MONTH_RAW")

echo
echo "=================================================="
echo "  Settings"
echo "=================================================="
echo "  Directory        : $RST_DIR"
echo "  Recurse          : $RECURSE"
echo "  Year range       : $START_YEAR – $END_YEAR"
echo "  Month to retain  : $RETAIN_MONTH"
echo "  Staging dir      : $(pwd)/restarts"
echo "  Dry run          : $DRY_RUN"
echo

# ================================================================
# Scan all .nc files in the directory
# ================================================================
echo "Scanning for restart files ..."
echo

FIND_DEPTH=(-maxdepth 1)
$RECURSE && FIND_DEPTH=()  # no depth limit if recursive

# Associative arrays:
#   comp_year_months[comp:year] = "01 03 07 ..."  (space-sep list of months seen)
declare -A comp_year_months

FILES_TO_REMOVE=()
FILES_TO_KEEP=()
TOTAL_BYTES=0

while IFS= read -r -d '' filepath; do
  fname=$(basename "$filepath")

  # Extract the first YYYY-MM-DD pattern from the filename.
  # This matches both CIME-style (YYYY-MM-DD-SSSSS) and
  # MPAS-style (YYYY-MM-DD_HH.MM.SS) dates.
  # The "|| true" prevents set -e from exiting when grep finds no match
  # (e.g. monthly history files with YYYY-MM-only dates are silently skipped).
  date_str=$(echo "$fname" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -1 || true)
  [[ -z "$date_str" ]] && continue

  file_year="${date_str:0:4}"
  file_month="${date_str:5:2}"

  # Skip files outside the requested year range
  (( file_year < START_YEAR || file_year > END_YEAR )) && continue

  comp=$(get_component "$fname")

  # Skip history files, initial-condition files, and anything else that
  # is not a restart (token after component must be r, rs, rst, or rh<N>)
  is_restart_file "$fname" "$comp" || continue

  key="${comp}:${file_year}"

  # Record which months exist for this component+year
  if [[ -n "${comp_year_months[$key]+x}" ]]; then
    comp_year_months[$key]+=" $file_month"
  else
    comp_year_months[$key]="$file_month"
  fi

  # Queue files for removal or retention
  if [[ "$file_month" != "$RETAIN_MONTH" ]]; then
    FILES_TO_REMOVE+=("$filepath")
    fsize=$(stat -c '%s' "$filepath" 2>/dev/null || echo 0)
    TOTAL_BYTES=$(( TOTAL_BYTES + fsize ))
  else
    FILES_TO_KEEP+=("$filepath")
  fi

done < <(find "$RST_DIR" "${FIND_DEPTH[@]}" -name "*.nc" -print0 | sort -z)

# ================================================================
# Month availability check: warn if retain month is missing for
# any (component, year) pair that has files in the year range.
# ================================================================
echo "--------------------------------------------------"
echo "  Month availability check (retain month = $RETAIN_MONTH)"
echo "--------------------------------------------------"
MISSING_FOUND=false
for key in $(printf '%s\n' "${!comp_year_months[@]}" | sort); do
  comp="${key%%:*}"
  year="${key##*:}"
  months="${comp_year_months[$key]}"
  if ! printf '%s\n' $months | sort -u | grep -qx "$RETAIN_MONTH"; then
    avail=$(printf '%s\n' $months | sort -u | tr '\n' ' ')
    echo "  WARNING: month $RETAIN_MONTH not found for component '$comp', year $year"
    echo "           Available months this year: ${avail}"
    MISSING_FOUND=true
  fi
done
$MISSING_FOUND || echo "  OK — month $RETAIN_MONTH is present for all component/year combinations found."
echo

# ================================================================
# Dry run listing  or  actual deletion
# ================================================================
N_FILES=${#FILES_TO_REMOVE[@]}

if [[ $N_FILES -eq 0 ]]; then
  echo "No files to move. Exiting."
  exit 0
fi

SIZE_FMT=$(format_bytes "$TOTAL_BYTES")

STAGE_DIR="$(pwd)/restarts"

if $DRY_RUN; then
  echo "--------------------------------------------------"
  echo "  DRY RUN — files that WOULD be moved to ./restarts/"
  echo "  Total: $N_FILES files, $SIZE_FMT"
  echo "--------------------------------------------------"
  for f in "${FILES_TO_REMOVE[@]}"; do
    fsize=$(stat -c '%s' "$f" 2>/dev/null || echo 0)
    printf "  %10s  %s\n" "$(format_bytes "$fsize")" "$(basename "$f")"
  done
  echo
  echo "  Summary: $N_FILES files, $SIZE_FMT would be moved to: $STAGE_DIR"
  echo
  N_KEEP=${#FILES_TO_KEEP[@]}
  echo "--------------------------------------------------"
  echo "  DRY RUN — files that WOULD be retained ($N_KEEP)"
  echo "--------------------------------------------------"
  for f in "${FILES_TO_KEEP[@]}"; do
    fsize=$(stat -c '%s' "$f" 2>/dev/null || echo 0)
    printf "  %10s  %s\n" "$(format_bytes "$fsize")" "$(basename "$f")"
  done
  echo
  echo "(Re-run and answer N to dry-run prompt to perform the actual move.)"
else
  echo "--------------------------------------------------"
  echo "  Ready to move $N_FILES files ($SIZE_FMT) to ./restarts/"
  echo "--------------------------------------------------"
  echo "  Destination: $STAGE_DIR"
  echo
  read -rp "  Confirm move? [y/N]: " CONFIRM
  if [[ ! "${CONFIRM,,}" =~ ^y ]]; then
    echo "Aborted. No files were moved."
    exit 0
  fi
  mkdir -p "$STAGE_DIR"
  MOVED=0
  TOTAL_MOVED=0
  for f in "${FILES_TO_REMOVE[@]}"; do
    fsize=$(stat -c '%s' "$f" 2>/dev/null || echo 0)
    echo "  Moving: $(basename "$f")"
    mv "$f" "$STAGE_DIR/"
    (( MOVED++ )) || true
    TOTAL_MOVED=$(( TOTAL_MOVED + fsize ))
  done
  echo
  echo "Done. Moved $MOVED files ($(format_bytes "$TOTAL_MOVED")) to: $STAGE_DIR"
fi
