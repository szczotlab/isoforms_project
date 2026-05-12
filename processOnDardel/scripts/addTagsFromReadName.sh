#!/usr/bin/env bash
set -euo pipefail

awk '
BEGIN {OFS="\t"}
# Keep header lines as-is
/^@/ {print; next}
{
	cb=""; ub="";
	# Extract CB and UB from read name
	if ($1 ~  /_OX:[A-Z]+_CR:[A-Z-]+$/) {
		match($1, /_OX:([A-Z]+)_CR:([A-Z-]+)$/, arr)
		ub=arr[1]; cb=arr[2] 
		$0 = $0 "\tCR:Z:" cb "\tOX:Z:" ub
	}
	else if ($1 ~ /_CR:[A-Z-]+$/) {
		match($1, /_CR:([A-Z-]+)$/, arr)
		cb=arr[1]
		$0 = $0 "\tCR:Z:" cb
	}
	print
}'

