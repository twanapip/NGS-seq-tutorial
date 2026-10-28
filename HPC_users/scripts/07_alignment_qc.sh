#!/usr/bin/env bash
# Alignment QC: flagstat, alignment summary, insert size, and coverage (WGS metrics or WES hybrid-selection metrics).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
M="$PROJ/qc/metrics"
for S in $(samples_from_args "$@"); do
  B="$PROJ/bam/$S.recal.bam"
  samtools flagstat "$B" > "$M/$S.flagstat.txt"
  gatkj CollectAlignmentSummaryMetrics -R "$REF_FASTA" -I "$B" -O "$M/$S.alignment_summary.txt"
  gatkj CollectInsertSizeMetrics -I "$B" -O "$M/$S.insert_size.txt" -H "$M/$S.insert_size_hist.pdf" \
    || echo "WARN: insert-size histogram failed (needs R); metrics table still written"
  if [ "$SEQ_TYPE" = "WES" ]; then
    gatkj CollectHsMetrics -R "$REF_FASTA" -I "$B" -O "$M/$S.hs_metrics.txt" \
      --BAIT_INTERVALS "$TARGETS_IL" --TARGET_INTERVALS "$TARGETS_IL" --PER_TARGET_COVERAGE "$M/$S.per_target_coverage.txt"
  else
    gatkj CollectWgsMetrics -R "$REF_FASTA" -I "$B" -O "$M/$S.wgs_metrics.txt"
  fi
done
