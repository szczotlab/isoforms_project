#!/usr/bin/env bash
set -euo pipefail

awk '
BEGIN { OFS="\t" }

# Keep header lines as-is
/^@/ { print; next }

{
    # Loop over all fields to find CR tag
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^CR:Z:/) {
            # Extract CR sequence
            match($i, /^CR:Z:(.*)$/, arr)
            cr = arr[1]

            # Replace N with A
            gsub(/N/, "A", cr)

            # Rebuild the field
            $i = "CR:Z:" cr
        }
    }

    print
}' 

