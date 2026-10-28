#!/usr/bin/env bash
# BWA-MEM alignment (read group included, GATK needs it) piped into a sorted, indexed BAM.
# Usage: bash scripts/05_align.sh [ID ...]
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env core
for S in $(samples_from_args "$@"); do
  if [ "$USE_TRIMMED" = 1 ]; then
    R1="$PROJ/data/trimmed/${S}_R1_val_1.fq.gz"; R2="$PROJ/data/trimmed/${S}_R2_val_2.fq.gz"
  else
    R1="$PROJ/data/raw/${S}_R1.fastq.gz"; R2="$PROJ/data/raw/${S}_R2.fastq.gz"
  fi
  [ -s "$R1" ] && [ -s "$R2" ] || { echo "Reads not found for $S: $R1"; exit 1; }
  bwa mem -t "$CORES" -M -R "@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA\tLB:${S}_lib" "$REF_FASTA" "$R1" "$R2" \
    | samtools sort -@ "$CORES" -T "$TMP_DIR/${S}" -o "$PROJ/bam/${S}.sorted.bam" -
  samtools index "$PROJ/bam/${S}.sorted.bam"
done
