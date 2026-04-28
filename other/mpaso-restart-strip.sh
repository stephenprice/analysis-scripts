#!/bin/bash

source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
set -euo pipefail

infile="mpaso-restart.nc"
outfile="mpaso-restart-mesh.nc"

[[ -f "$infile" ]] || { echo "ERROR: $infile not found"; exit 1; }

time_vars=$(
  ncdump -h "$infile" | awk '
    /^[[:space:]]*(byte|char|short|int|float|double|ubyte|ushort|uint|int64|uint64|string)[[:space:]]+[A-Za-z_][A-Za-z0-9_]*\(/ {
      line=$0
      sub(/^[[:space:]]*(byte|char|short|int|float|double|ubyte|ushort|uint|int64|uint64|string)[[:space:]]+/, "", line)
      name=line
      sub(/\(.*/, "", name)

      dims=line
      sub(/^[^(]*\(/, "", dims)
      sub(/\).*/, "", dims)
      gsub(/[[:space:]]/, "", dims)

      n=split(dims, a, ",")
      for (i=1; i<=n; i++) {
        if (a[i] == "Time") {
          print name
          break
        }
      }
    }
  ' | paste -sd, -
)

[[ -n "${time_vars}" ]] || { echo "ERROR: no variables with dimension Time were found"; exit 1; }

ncks -O -x -v "$time_vars" "$infile" "$outfile"
