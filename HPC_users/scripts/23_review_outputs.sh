#!/usr/bin/env bash
# Build the MultiQC report and show what each stage produced.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env core
multiqc -f -o "$PROJ/qc/multiqc" "$PROJ/qc" "$PROJ/data/trimmed" 2>/dev/null || echo "MultiQC found nothing yet"
echo "== QC summaries"; for f in "$PROJ"/qc/*.qc_summary.txt; do [ -e "$f" ] && { cat "$f"; echo; }; done
for d in bam gvcf vcf cohort annotation somatic cnv maf gwas; do
  n=$(find "$PROJ/$d" -type f 2>/dev/null | wc -l); echo "$d: $n files"
done
echo "MultiQC report: $PROJ/qc/multiqc/multiqc_report.html"
