#!/usr/bin/env bash
# Mark duplicates, then Base Quality Score Recalibration -> bam/<ID>.recal.bam
# WES: BQSR statistics are learned on the targets only, but the final BAM keeps ALL reads (CNV needs off-target reads).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
for S in $(samples_from_args "$@"); do
  B="$PROJ/bam/$S"
  gatkj MarkDuplicates -I "$B.sorted.bam" -O "$B.dedup.bam" -M "$PROJ/qc/metrics/$S.dup_metrics.txt" \
    --CREATE_INDEX true --TMP_DIR "$TMP_DIR"
  gatkj BaseRecalibrator -R "$REF_FASTA" -I "$B.dedup.bam" --known-sites "$DBSNP" --known-sites "$MILLS" \
    "${INTERVAL_ARGS[@]}" -O "$B.recal.table"
  gatkj ApplyBQSR -R "$REF_FASTA" -I "$B.dedup.bam" --bqsr-recal-file "$B.recal.table" -O "$B.recal.bam"
  if [ "$CLEAN_INTERMEDIATE" = 1 ]; then rm -f "$B.sorted.bam" "$B.sorted.bam.bai" "$B.dedup.bam" "$B.dedup.bai" "$B.dedup.bam.bai"; fi
done
