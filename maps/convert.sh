#!/bin/bash

# Check if the input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file="$1"

# Function to convert degrees, minutes, seconds to decimal degrees
convert_to_decimal() {
    local degrees=$1
    local minutes=$2
    local seconds=$3
    echo "scale=6; $degrees + ($minutes / 60) + ($seconds / 3600)" | bc
}

# Process the file line by line
while read -r lat_deg lat_min lat_sec lon_deg lon_min lon_sec; do
    # Convert latitude and longitude to decimal degrees
    decimal_lat=$(convert_to_decimal "$lat_deg" "$lat_min" "$lat_sec")
    decimal_lon=$(convert_to_decimal "$lon_deg" "$lon_min" "$lon_sec")

    # Print the result
    echo "$decimal_lat, $decimal_lon"
done < "$input_file"
