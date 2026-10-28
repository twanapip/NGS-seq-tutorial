#!/usr/bin/env bash
# Create the project folders and empty list files you need to fill in.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
mkdir -p "$PROJ"/{lists,data/raw,data/trimmed,qc/fastqc_raw,qc/metrics,bam,gvcf,vcf,cohort,annotation,somatic,cnv,maf,gwas,logs,tmp,targets}
mk() { [ -e "$1" ] || printf '%s\n' "$2" > "$1"; }
mk "$SAMPLE_LIST"   "# One sample ID per line. FASTQs must be data/raw/<ID>_R1.fastq.gz and data/raw/<ID>_R2.fastq.gz"
mk "$PAIRS_FILE"    "# tumor_id<TAB>normal_id   (use - as normal_id for tumor-only). Tumor IDs must be unique."
mk "$NORMALS_FILE"  "# Normal IDs for pooled normal / custom panel of normals / CNV reference"
mk "$CASES_FILE"    "# Patient (case) IDs for disease/patient-specific variants"
mk "$CONTROLS_FILE" "# Control IDs for disease/patient-specific variants"
echo "Project: $PROJ"
echo "1) Put FASTQs in $PROJ/data/raw   2) Fill the files in $PROJ/lists"
