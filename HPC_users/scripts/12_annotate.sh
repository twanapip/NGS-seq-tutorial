#!/usr/bin/env bash
# Annotate one VCF with every tool in ANNOTATORS.
# Usage: bash scripts/12_annotate.sh <sample_id | path/to.vcf.gz> [germline|somatic]
#   germline + ID  -> vcf/<ID>.pass.vcf.gz          somatic + tumor ID -> somatic/<MODE>/<ID>.high_confidence.vcf.gz
# Output: annotation/<vcf name>/{snpeff,vep,funcotator,annovar}*
set -euo pipefail
source "$(dirname "$0")/../config.sh"
A="${1:?usage: 12_annotate.sh <id|vcf.gz> [germline|somatic]}"; KIND="${2:-germline}"
if [ ! -f "$A" ]; then
  if [ "$KIND" = somatic ]; then A="$PROJ/somatic/$SOMATIC_MODE/$A.high_confidence.vcf.gz"; else A="$PROJ/vcf/$A.pass.vcf.gz"; fi
fi
[ -s "$A" ] || { echo "VCF not found: $A"; exit 1; }
NAME="$(basename "$A" .vcf.gz)"; OUT="$PROJ/annotation/$NAME"; mkdir -p "$OUT"

if has_annotator snpeff; then
  use_env snpeff
  snpEff -Xmx"${SNPEFF_MEM}"g ann -v -canon -hgvs -stats "$OUT/snpeff_stats.html" -csvStats "$OUT/snpeff_stats.csv" "$SNPEFF_DB" "$A" \
    | bgzip -c > "$OUT/snpeff.vcf.gz"
  tabix -f -p vcf "$OUT/snpeff.vcf.gz"
fi

if has_annotator funcotator; then
  use_env gatk
  if [ "$KIND" = somatic ]; then DS="$(funco_dir s)"; else DS="$(funco_dir g)"; fi
  [ -n "$DS" ] || { echo "Funcotator data sources missing: run scripts/11_annotation_db.sh"; exit 1; }
  for FMT in VCF MAF; do
    EXT=$([ "$FMT" = VCF ] && echo vcf || echo maf)
    gatk Funcotator --variant "$A" --reference "$REF_FASTA" --ref-version "$FUNCO_REF" --data-sources-path "$DS" \
      --output "$OUT/funcotator.$EXT" --output-file-format "$FMT" --remove-filtered-variants false --transcript-selection-mode CANONICAL
  done
fi

if has_annotator vep; then
  use_env vep
  VA=(--input_file "$A" --format vcf --vcf --output_file "$OUT/vep.vcf" --force_overwrite
      --cache --offline --dir_cache "$VEP_CACHE" --cache_version "$VEP_CACHE_VERSION" --assembly "$VEP_ASSEMBLY" --species homo_sapiens
      --fasta "$REF_FASTA" --symbol --canonical --biotype --hgvs --hgvsg --protein --ccds --uniprot --domains --numbers --tsl
      --regulatory --terms SO --no_escape --af --af_1kg --af_gnomade --af_gnomadg --max_af
      --fork "$CORES" --stats_file "$OUT/vep_stats.html" --warning_file "$OUT/vep_warnings.txt")
  if [ -n "$DBNSFP_FILE" ]; then
    VA+=(--dir_plugins "$VEP_CACHE/Plugins" --plugin "dbNSFP,$DBNSFP_FILE,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_pred,CADD_raw,CADD_phred")
  fi
  vep "${VA[@]}"
  bgzip -f "$OUT/vep.vcf"; tabix -f -p vcf "$OUT/vep.vcf.gz"
fi

if has_annotator annovar; then
  if [ -z "$ANNOVAR_DIR" ] || [ ! -f "$ANNOVAR_DIR/table_annovar.pl" ]; then
    echo "WARN: ANNOVAR_DIR not set - skipping ANNOVAR"
  else
    use_env gatk   # only for bcftools
    mkdir -p "$OUT/annovar"
    bcftools view "$A" > "$OUT/annovar/input.vcf"
    perl "$ANNOVAR_DIR/convert2annovar.pl" -format vcf4 "$OUT/annovar/input.vcf" > "$OUT/annovar/input.avinput"
    perl "$ANNOVAR_DIR/table_annovar.pl" "$OUT/annovar/input.avinput" "$ANNOVAR_DIR/humandb/" -buildver "$ANNOVAR_BUILD" \
      -out "$OUT/annovar/annovar" -remove -protocol "$ANNOVAR_PROTOCOLS" -operation "$ANNOVAR_OPERATIONS" -nastring . -csvout
    rm -f "$OUT/annovar/input.vcf"
  fi
fi
echo "Annotated: $OUT"
