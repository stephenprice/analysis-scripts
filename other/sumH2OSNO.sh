#!/bin/bash

set -euo pipefail

# This file uses NCO to isolate and sum up the liquid water equivalent in the
# snowpack over Greenland from ELM history files using a Greenland mask.

mask='/global/cfs/cdirs/e3sm/sprice/temp/masks/glcmaskdata_0.5x0.5_GIS.20260310.nc'

usage() {
	echo "Usage: $0 [start_year end_year]"
	echo "If years are not provided, the script prompts for them interactively."
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	usage
	exit 0
fi

# Auto-discover case from any file matching *.elm.h0.*.nc in this directory.
first_file=$(ls -1 *.elm.h0.*.nc 2>/dev/null | head -n 1 || true)
if [[ -z "${first_file}" ]]; then
	echo "ERROR: No files matching '*.elm.h0.*.nc' found in $(pwd)"
	exit 1
fi
case=${first_file%%.elm.h0.*}

if [[ $# -ge 2 ]]; then
	start_year="$1"
	end_year="$2"
else
	read -r -p "Enter start year (e.g., 2012): " start_year
	read -r -p "Enter end year   (e.g., 2015): " end_year
fi

if ! [[ "${start_year}" =~ ^[0-9]{4}$ && "${end_year}" =~ ^[0-9]{4}$ ]]; then
	echo "ERROR: start_year and end_year must both be 4-digit years."
	exit 1
fi

if (( end_year < start_year )); then
	echo "ERROR: end_year (${end_year}) must be >= start_year (${start_year})."
	exit 1
fi

fileout="h2osno_${start_year}-${end_year}.nc"

inputs=()
for ((yy=start_year; yy<=end_year; yy++)); do
	while IFS= read -r f; do
		inputs+=("${f}")
	done < <(ls -1 "${case}.elm.h0.${yy}"*.nc 2>/dev/null || true)
done

if (( ${#inputs[@]} == 0 )); then
	echo "ERROR: No files found for case='${case}' and years ${start_year}-${end_year}."
	exit 1
fi

echo "Discovered case: ${case}"
echo "Input year range: ${start_year}-${end_year}"
echo "Number of files to concatenate: ${#inputs[@]}"
echo "Output file: ${fileout}"

# Source module stack that provides nco and python.
set +u
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
set -u

# Concatenate selected files while extracting required fields.
ncrcat -v H2OSNO,area "${inputs[@]}" "${fileout}"

# Copy the mask file and update dim names to match ELM history files.
cp "${mask}" ./gis-mask.nc
ncrename -d lsmlat,lat gis-mask.nc
ncrename -d lsmlon,lon gis-mask.nc

# Append the mask field that isolates Greenland.
ncks -A -v GLCMASK gis-mask.nc "${fileout}"
rm gis-mask.nc

# Convert H2OSNO to km, multiply by area and mask -> km^3 (Gt).
ncap2 -A -s 'maskedH2OSNOkm=H2OSNO/1000/1000*area*GLCMASK' "${fileout}"
ncatted -O -a long_name,maskedH2OSNOkm,o,c,"snow depth (liquid water) multipled by cell area and mask" "${fileout}"
ncatted -O -a units,maskedH2OSNOkm,o,c,"cubic km or Gt" "${fileout}"

# Area-integrated time series.
ncap2 -A -s 'sumMaskedH2OSNOkm=maskedH2OSNOkm.total($lat,$lon);' "${fileout}"
ncatted -O -a long_name,sumMaskedH2OSNOkm,o,c,"snow depth (liquid water) multipled by mask and cell area, summed over all cells" "${fileout}"
ncatted -O -a units,sumMaskedH2OSNOkm,o,c,"cubic km or Gt" "${fileout}"

# Convert total Gt to approximate global mean sea-level equivalent (mm).
ncap2 -A -s 'mmSeaLevelEquiv=sumMaskedH2OSNOkm/360' "${fileout}"
ncatted -O -a long_name,mmSeaLevelEquiv,o,c,"snow depth (liquid water) multipled by mask and cell area, summed over all cells and converted to approx. global sea-level equiv." "${fileout}"
ncatted -O -a units,mmSeaLevelEquiv,o,c,"mm of global sea-level equivalent" "${fileout}"

# Sea-level-equivalent change relative to first time sample.
ncap2 -A -s 'mmSeaLevelEquivChange=mmSeaLevelEquiv-mmSeaLevelEquiv(0)' "${fileout}"
ncatted -O -a long_name,mmSeaLevelEquivChange,o,c,"change in global sea-level equivalent relative to first time sample" "${fileout}"
ncatted -O -a units,mmSeaLevelEquivChange,o,c,"mm of global sea-level equivalent" "${fileout}"

echo "Done. Wrote ${fileout}"
