#!/usr/bin/env bash
# Submit a stage to SLURM with dependencies. Usage: ./submit_pipeline.sh <stage>
#   germline | cohort | somatic | cnv | disease | maf | gwas | review
# One-time setup first, on a login/transfer node (internet needed):
#   bash scripts/00_conda_envs.sh; bash scripts/01_reference.sh; bash scripts/03_project_dirs.sh
#   WES: bash scripts/02_targets.sh        annotators: bash scripts/11_annotation_db.sh
# Then put FASTQs in $PROJ/data/raw and fill $PROJ/lists/*.
set -euo pipefail
cd "$(dirname "$0")"
source ./config.sh
mkdir -p "$PROJ/logs"
N=$(read_list "$SAMPLE_LIST" | wc -l); NP=$(read_list "$PAIRS_FILE" | wc -l)
MEM=$(( MEM_GB + 8 ))

# J <name> <wrapper> <cpus> <array-size|""> <dependency|""> <script> [args...]  -> prints job id
J() {
  local name=$1 wrap=$2 cpus=$3 arr=$4 dep=$5 script=$6; shift 6
  local log="$PROJ/logs/${name}_%j.out"; [ -n "$arr" ] && log="$PROJ/logs/${name}_%A_%a.out"
  local o=(--parsable --job-name="$name" --partition="$PARTITION" --time="$TIME" --cpus-per-task="$cpus" --mem="${MEM}G"
           --output="$log" --mail-user="$MAIL" --mail-type=FAIL)
  [ -n "$arr" ] && o+=(--array="1-$arr")
  [ -n "$dep" ] && o+=(--dependency="$dep")
  sbatch "${o[@]}" --export=ALL,STEP="scripts/$script" "scripts/$wrap" "$@"
}
need() { [ "$1" -ge 1 ] || { echo "$2 is empty - fill it in first"; exit 1; }; }

STAGE="${1:-}"
case "$STAGE" in
  germline)
    need "$N" "samples.txt"
    a=$(J fastqc_trim slurm_array.sh "$CORES" "$N" "" 04_fastqc_trim.sh)
    b=$(J align       slurm_array.sh "$CORES" "$N" "aftercorr:$a" 05_align.sh)
    c=$(J markdup_bqsr slurm_array.sh "$CORES" "$N" "aftercorr:$b" 06_markdup_bqsr.sh)
    d=$(J align_qc    slurm_array.sh "$CPUS_LIGHT" "$N" "aftercorr:$c" 07_alignment_qc.sh)
    e=$(J call_germline slurm_array.sh "$CORES" "$N" "aftercorr:$c" 08_call_germline.sh)
    f=$(J filter      slurm_array.sh "$CPUS_LIGHT" "$N" "aftercorr:$e" 09_filter_variants.sh)
    g=$(J annotate    slurm_array.sh "$CORES" "$N" "aftercorr:$f" 12_annotate.sh germline)
    h=$(J qc_summary  slurm_array.sh "$CPUS_LIGHT" "$N" "aftercorr:$f,aftercorr:$d" 13_variant_stats_qc.sh)
    echo "Submitted germline chain; last jobs: annotate=$g qc=$h" ;;
  cohort)
    need "$N" "samples.txt"
    a=$(J joint slurm_single.sh "$CORES" "" "" 10_joint_calling.sh)
    b=$(J annotate_cohort slurm_single.sh "$CORES" "" "afterok:$a" 12_annotate.sh "$PROJ/cohort/cohort.pass.vcf.gz" germline)
    echo "Submitted: joint=$a annotate=$b" ;;
  somatic)
    need "$NP" "pairs.tsv"
    a=$(J somatic_prep slurm_single.sh "$CORES" "" "" 14_somatic_prepare.sh)
    b=$(J mutect2 slurm_pairs.sh "$CORES" "$NP" "afterok:$a" 15_somatic_mutect2.sh)
    c=$(J somatic_filter slurm_pairs.sh "$CPUS_LIGHT" "$NP" "aftercorr:$b" 16_somatic_filter.sh)
    d=$(J somatic_annot slurm_pairs.sh "$CORES" "$NP" "aftercorr:$c" 12_annotate.sh somatic)
    echo "Submitted: prep=$a mutect2=$b filter=$c annotate=$d" ;;
  cnv)     J cnv slurm_single.sh "$CORES" "" "" 17_cnv.sh "$CNV_MODE" ;;
  disease) J disease slurm_single.sh "$CPUS_LIGHT" "" "" 18_disease_variants.sh ;;
  maf)     J maf slurm_single.sh "$CPUS_LIGHT" "" "" 19_maf_visualize.sh ;;
  gwas)
    a=$(J gwas_prepare slurm_single.sh "$CORES" "" "" 20_gwas_prepare.sh)
    b=$(J gwas_qc_pca slurm_single.sh "$CORES" "" "afterok:$a" 21_gwas_qc_pca.sh)
    c=$(J gwas_assoc slurm_single.sh "$CORES" "" "afterok:$b" 22_gwas_assoc_plot.sh)
    echo "Submitted: prepare=$a qc_pca=$b assoc=$c" ;;
  review)  J review slurm_single.sh "$CPUS_LIGHT" "" "" 23_review_outputs.sh ;;
  *) sed -n '2,6p' "$0"; exit 1 ;;
esac
