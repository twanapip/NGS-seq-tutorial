#!/usr/bin/env bash
# FastQC on raw reads, Trim Galore (adapters, bases below Q20, reads under 50 bp), FastQC on trimmed reads.
# Usage: bash scripts/04_fastqc_trim.sh [ID ...]     (no argument = every ID in samples.txt)
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env core
TC=$(( CORES < 4 ? CORES : 4 ))
for S in $(samples_from_args "$@"); do
  R1="$PROJ/data/raw/${S}_R1.fastq.gz"; R2="$PROJ/data/raw/${S}_R2.fastq.gz"
  [ -s "$R1" ] && [ -s "$R2" ] || { echo "Missing FASTQ for $S in $PROJ/data/raw/"; exit 1; }
  fastqc -t "$CORES" -o "$PROJ/qc/fastqc_raw" "$R1" "$R2"
  trim_galore --paired --quality 20 --length 50 --cores "$TC" --fastqc --output_dir "$PROJ/data/trimmed" "$R1" "$R2"
done
