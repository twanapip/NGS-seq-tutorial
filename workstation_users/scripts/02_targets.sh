#!/usr/bin/env bash
# WES only: check the capture-kit BED matches the genome and build the GATK interval list.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
[ "$SEQ_TYPE" = "WES" ] || { echo "SEQ_TYPE is not WES - nothing to do"; exit 0; }
use_env gatk
[ -s "$TARGETS" ] || { echo "TARGETS file not found: $TARGETS"; exit 1; }
BAD=$(awk 'NR==FNR{c[$1]=1;next} !/^(track|browser|#)/ && !($1 in c){print $1}' "$REF_FASTA.fai" "$TARGETS" | sort -u | head -5)
if [ -n "$BAD" ]; then
  echo "Contigs in TARGETS missing from the reference (chr1 vs 1 mismatch or wrong build?):"; echo "$BAD"; exit 1
fi
mkdir -p "$PROJ/targets"
gatk BedToIntervalList -I "$TARGETS" -O "$TARGETS_IL" -SD "$REF_DICT"
echo "Interval list: $TARGETS_IL"
