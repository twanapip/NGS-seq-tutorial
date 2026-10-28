#!/usr/bin/env bash
# One-time prep for the chosen SOMATIC_MODE (config.sh):
#   pooled     -> merge the normals in normals.txt into one BAM with a single read group "pooled_normal"
#   custom_pon -> build a Panel of Normals from normals.txt (use 40+ normals from the same kit/platform in real studies)
#   matched / tumor_only -> nothing to prepare
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
mapfile -t NORMALS < <(read_list "$NORMALS_FILE")
case "$SOMATIC_MODE" in
  pooled)
    [ "${#NORMALS[@]}" -ge 2 ] || { echo "pooled mode needs 2+ IDs in normals.txt"; exit 1; }
    D="$PROJ/somatic/pooled"; mkdir -p "$D"
    BAMS=(); for N in "${NORMALS[@]}"; do BAMS+=("$PROJ/bam/$N.recal.bam"); done
    samtools merge -f -@ "$CORES" "$D/merged.bam" "${BAMS[@]}"
    # Mutect2 looks the normal up by read-group sample name, so the merged BAM must carry ONE sample name
    gatkj AddOrReplaceReadGroups -I "$D/merged.bam" -O "$D/pooled_normal.bam" --RGID pooled_normal --RGSM pooled_normal \
      --RGLB pooled --RGPL ILLUMINA --RGPU pooled --CREATE_INDEX true
    rm -f "$D/merged.bam" "$D/merged.bam.bai"
    echo "Pooled normal: $D/pooled_normal.bam" ;;
  custom_pon)
    [ "${#NORMALS[@]}" -ge 2 ] || { echo "custom_pon needs normals in normals.txt"; exit 1; }
    [ "${#NORMALS[@]}" -ge 40 ] || echo "WARN: ${#NORMALS[@]} normals - a useful PoN needs ~40+. Fine for learning, weak for real work."
    D="$PROJ/somatic/custom_pon"; mkdir -p "$D/normals_calls"
    for N in "${NORMALS[@]}"; do
      gatkj Mutect2 -R "$REF_FASTA" -I "$PROJ/bam/$N.recal.bam" --max-mnp-distance 0 "${INTERVAL_ARGS[@]}" \
        -O "$D/normals_calls/$N.vcf.gz"
    done
    PARTS=()
    for U in "${UNITS[@]}"; do
      set_unit "$U"; VARGS=(); for N in "${NORMALS[@]}"; do VARGS+=(-V "$D/normals_calls/$N.vcf.gz"); done
      rm -rf "$D/db_$U"
      gatkj GenomicsDBImport -R "$REF_FASTA" --genomicsdb-workspace-path "$D/db_$U" --tmp-dir "$TMP_DIR" "${UNIT_L[@]}" "${VARGS[@]}"
      gatkj CreateSomaticPanelOfNormals -R "$REF_FASTA" -V "gendb://$D/db_$U" --germline-resource "$GNOMAD" -O "$D/pon_$U.vcf.gz"
      rm -rf "$D/db_$U"; PARTS+=(-I "$D/pon_$U.vcf.gz")
    done
    gatkj MergeVcfs "${PARTS[@]}" -O "$D/custom_pon.vcf.gz"
    rm -f "$D"/pon_*.vcf.gz*
    echo "Custom PoN: $D/custom_pon.vcf.gz" ;;
  matched|tumor_only) echo "SOMATIC_MODE=$SOMATIC_MODE needs no preparation" ;;
  *) echo "SOMATIC_MODE must be matched | pooled | custom_pon | tumor_only"; exit 1 ;;
esac
