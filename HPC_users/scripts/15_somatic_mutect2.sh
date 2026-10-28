#!/usr/bin/env bash
# Mutect2 for each tumor in pairs.tsv, using SOMATIC_MODE. Output: somatic/<MODE>/<tumor>.raw.vcf.gz (+ .stats, f1r2)
# Usage: bash scripts/15_somatic_mutect2.sh [tumor_id ...]
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
OUT="$PROJ/somatic/$SOMATIC_MODE"; mkdir -p "$OUT"
for T in $(tumors_from_args "$@"); do
  N="$(normal_of "$T")"
  ARGS=(-R "$REF_FASTA" -I "$PROJ/bam/$T.recal.bam" -tumor "$T" --germline-resource "$GNOMAD"
        --f1r2-tar-gz "$OUT/$T.f1r2.tar.gz" -O "$OUT/$T.raw.vcf.gz"
        --native-pair-hmm-threads "$HMM_THREADS" --max-reads-per-alignment-start 50 "${INTERVAL_ARGS[@]}")
  case "$SOMATIC_MODE" in
    matched)
      { [ -n "$N" ] && [ "$N" != "-" ]; } || { echo "matched mode: no normal for $T in pairs.tsv"; exit 1; }
      ARGS+=(-I "$PROJ/bam/$N.recal.bam" -normal "$N" --panel-of-normals "$PON") ;;
    pooled)
      [ -s "$PROJ/somatic/pooled/pooled_normal.bam" ] || { echo "Run scripts/14_somatic_prepare.sh first"; exit 1; }
      ARGS+=(-I "$PROJ/somatic/pooled/pooled_normal.bam" -normal pooled_normal --panel-of-normals "$PON") ;;
    custom_pon)
      [ -s "$PROJ/somatic/custom_pon/custom_pon.vcf.gz" ] || { echo "Run scripts/14_somatic_prepare.sh first"; exit 1; }
      ARGS+=(--panel-of-normals "$PROJ/somatic/custom_pon/custom_pon.vcf.gz" --max-mnp-distance 0) ;;
    tumor_only)
      ARGS+=(--panel-of-normals "$PON" --max-mnp-distance 0) ;;
    *) echo "SOMATIC_MODE must be matched | pooled | custom_pon | tumor_only"; exit 1 ;;
  esac
  gatkj Mutect2 "${ARGS[@]}"
done
