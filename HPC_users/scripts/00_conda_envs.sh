#!/usr/bin/env bash
# One conda env per tool group (separate envs avoid Java/R version conflicts, e.g. GATK vs SnpEff).
# Usage: bash scripts/00_conda_envs.sh [core gatk snpeff vep cnv gwas maf]     (no argument = all)
set -euo pipefail
source "$(dirname "$0")/../config.sh"
CH=(-c conda-forge -c bioconda)
mk() {
  local name="${ENV_PREFIX}_$1"; shift
  if conda env list | awk '{print $1}' | grep -qx "$name"; then echo "exists: $name"; else conda create -y -n "$name" "${CH[@]}" "$@"; fi
}
WANT=("$@"); [ "${#WANT[@]}" -gt 0 ] || WANT=(core gatk snpeff vep cnv gwas maf)
for e in "${WANT[@]}"; do
  case "$e" in
    core)   mk core fastqc multiqc trim-galore cutadapt bwa samtools bcftools htslib bedtools wget ;;
    gatk)   mk gatk gatk4 samtools bcftools htslib bedtools ;;
    snpeff) mk snpeff snpeff snpsift htslib ;;
    vep)    mk vep ensembl-vep samtools htslib git ;;
    cnv)    mk cnv cnvkit bioconductor-dnacopy bedtools bcftools htslib ;;
    gwas)   mk gwas plink gcta bcftools htslib r-base r-ggplot2 r-dplyr r-data.table ;;
    maf)    mk maf r-base r-tidyverse bioconductor-maftools bioconductor-mutationalpatterns \
               bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-bsgenome.hsapiens.ucsc.hg19 r-pheatmap r-rcolorbrewer ;;
    *) echo "unknown env: $e"; exit 1 ;;
  esac
done
echo "Done. Activate one with: conda activate ${ENV_PREFIX}_gatk"
