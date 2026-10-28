#!/usr/bin/env bash
# Joint genotyping of all GVCFs in samples.txt -> cohort/cohort.pass.vcf.gz
# WGS: one chromosome at a time (GenomicsDBImport); WES: the target file in one go.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
C="$PROJ/cohort"; mkdir -p "$C"
MAP="$C/sample_map.txt"; : > "$MAP"
for S in $(read_list "$SAMPLE_LIST"); do printf '%s\t%s\n' "$S" "$PROJ/gvcf/$S.g.vcf.gz" >> "$MAP"; done
[ "$(wc -l < "$MAP")" -ge 2 ] || { echo "Joint calling needs at least 2 samples"; exit 1; }
PARTS=()
for U in "${UNITS[@]}"; do
  set_unit "$U"; MERGE=(); [ "$U" = all ] && MERGE=(--merge-input-intervals)
  rm -rf "$C/db_$U"
  gatkj GenomicsDBImport --sample-name-map "$MAP" --genomicsdb-workspace-path "$C/db_$U" \
    --reader-threads "$HMM_THREADS" --tmp-dir "$TMP_DIR" "${UNIT_L[@]}" "${MERGE[@]}"
  gatkj GenotypeGVCFs -R "$REF_FASTA" -V "gendb://$C/db_$U" -O "$C/raw_$U.vcf.gz" "${UNIT_L[@]}"
  rm -rf "$C/db_$U"; PARTS+=(-I "$C/raw_$U.vcf.gz")
done
gatkj MergeVcfs "${PARTS[@]}" -O "$C/cohort.raw.vcf.gz"
rm -f "$C"/raw_*.vcf.gz*
filter_variants "$C/cohort.raw.vcf.gz" "$C/cohort"
echo "Cohort VCF: $C/cohort.pass.vcf.gz"
