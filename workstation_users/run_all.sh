#!/usr/bin/env bash
# Run a stage of the analysis in order. Usage: ./run_all.sh <stage>
#   setup     envs, reference, (WES targets), project folders, annotation databases   [internet; run once]
#   germline  FASTQ -> trimmed -> BAM -> recalibrated BAM -> QC -> variants -> filter -> annotate -> QC summary
#   cohort    joint calling of all samples + annotate the cohort VCF
#   somatic   Mutect2 (SOMATIC_MODE) -> filters -> annotate
#   cnv       CNVkit (CNV_MODE)
#   disease   patient-specific / disease-specific variants (run cohort first)
#   maf       somatic VCF -> MAF -> maftools plots
#   gwas      PLINK QC, PCA, association, plots
#   review    MultiQC + summary
# Before "germline": put FASTQs in $PROJ/data/raw and fill $PROJ/lists/samples.txt (see README).
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
run() { echo; echo "######## $*"; bash "scripts/$@"; }
STAGE="${1:-}"
case "$STAGE" in
  setup)
    run 00_conda_envs.sh; run 01_reference.sh; run 03_project_dirs.sh
    [ "$SEQ_TYPE" = WES ] && run 02_targets.sh
    run 11_annotation_db.sh ;;
  germline)
    run 04_fastqc_trim.sh; run 05_align.sh; run 06_markdup_bqsr.sh; run 07_alignment_qc.sh
    run 08_call_germline.sh; run 09_filter_variants.sh
    for S in $(read_list "$SAMPLE_LIST"); do run 12_annotate.sh "$S" germline; done
    run 13_variant_stats_qc.sh ;;
  cohort)
    run 10_joint_calling.sh; run 12_annotate.sh "$PROJ/cohort/cohort.pass.vcf.gz" germline ;;
  somatic)
    run 14_somatic_prepare.sh; run 15_somatic_mutect2.sh; run 16_somatic_filter.sh
    for T in $(tumors); do run 12_annotate.sh "$T" somatic; done ;;
  cnv)      run 17_cnv.sh "$CNV_MODE" ;;
  disease)  run 18_disease_variants.sh ;;
  maf)      run 19_maf_visualize.sh ;;
  gwas)     run 20_gwas_prepare.sh; run 21_gwas_qc_pca.sh; run 22_gwas_assoc_plot.sh ;;
  review)   run 23_review_outputs.sh ;;
  *) sed -n '2,15p' "$0"; exit 1 ;;
esac
echo; echo "Stage '$STAGE' finished."
