#!/usr/bin/env bash
# HaplotypeCaller in GVCF mode (keeps reference confidence so samples can be joint-called later), then genotype one sample.
# Scatter tip: for big WGS you can run with -L per chromosome in parallel and merge with MergeVcfs.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
for S in $(samples_from_args "$@"); do
  gatkj HaplotypeCaller -R "$REF_FASTA" -I "$PROJ/bam/$S.recal.bam" -O "$PROJ/gvcf/$S.g.vcf.gz" \
    -ERC GVCF --dbsnp "$DBSNP" --native-pair-hmm-threads "$HMM_THREADS" "${INTERVAL_ARGS[@]}"
  gatkj GenotypeGVCFs -R "$REF_FASTA" -V "$PROJ/gvcf/$S.g.vcf.gz" -O "$PROJ/vcf/$S.raw.vcf.gz" "${INTERVAL_ARGS[@]}"
done
