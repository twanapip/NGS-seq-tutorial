#!/usr/bin/env bash
# Contamination, read-orientation model, FilterMutectCalls, PASS extraction and mode-specific high-confidence filters.
# Output: somatic/<MODE>/<tumor>.{filtered,pass,high_confidence}.vcf.gz and <tumor>.mutations.tsv
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
OUT="$PROJ/somatic/$SOMATIC_MODE"; C="$OUT/contamination"; mkdir -p "$C"
PL=(-L "$COMMON"); [ -n "$TARGETS" ] && PL+=(-L "$TARGETS" --interval-set-rule INTERSECTION)
for T in $(tumors_from_args "$@"); do
  N="$(normal_of "$T")"; NB=""; NID=""
  case "$SOMATIC_MODE" in
    matched) NB="$PROJ/bam/$N.recal.bam"; NID="$N" ;;
    pooled)  NB="$PROJ/somatic/pooled/pooled_normal.bam"; NID="pooled_normal" ;;
  esac
  gatkj GetPileupSummaries -I "$PROJ/bam/$T.recal.bam" -V "$COMMON" "${PL[@]}" -O "$C/$T.pileups.table"
  FARGS=()
  if [ -n "$NB" ]; then
    gatkj GetPileupSummaries -I "$NB" -V "$COMMON" "${PL[@]}" -O "$C/$T.normal_pileups.table"
    gatkj CalculateContamination -I "$C/$T.pileups.table" -matched "$C/$T.normal_pileups.table" \
      -O "$C/$T.contamination.table" --tumor-segmentation "$C/$T.segments.table"
    FARGS+=(--tumor-segmentation "$C/$T.segments.table")
  else
    gatkj CalculateContamination -I "$C/$T.pileups.table" -O "$C/$T.contamination.table"
  fi
  gatkj LearnReadOrientationModel -I "$OUT/$T.f1r2.tar.gz" -O "$OUT/$T.orientation_model.tar.gz"
  gatkj FilterMutectCalls -R "$REF_FASTA" -V "$OUT/$T.raw.vcf.gz" --contamination-table "$C/$T.contamination.table" \
    --ob-priors "$OUT/$T.orientation_model.tar.gz" "${FARGS[@]}" -O "$OUT/$T.filtered.vcf.gz"
  bcftools view -f PASS -Oz -o "$OUT/$T.pass.vcf.gz" "$OUT/$T.filtered.vcf.gz"; tabix -f -p vcf "$OUT/$T.pass.vcf.gz"

  # sample positions in the VCF are looked up by name (never assume tumor=0)
  mapfile -t SN < <(bcftools query -l "$OUT/$T.pass.vcf.gz")
  TI=0; NI=""
  for i in "${!SN[@]}"; do [ "${SN[$i]}" = "$T" ] && TI=$i; { [ -n "$NID" ] && [ "${SN[$i]}" = "$NID" ]; } && NI=$i; done
  case "$SOMATIC_MODE" in
    matched)    AF=0.05; DP=10; TLOD=6.3 ;;
    pooled)     AF=0.06; DP=12; TLOD=6.3 ;;
    custom_pon) AF=0.05; DP=10; TLOD=6.3 ;;
    tumor_only) AF=0.10; DP=20; TLOD=10.0 ;;
  esac
  EXPR="FORMAT/AF[$TI:0] >= $AF && FORMAT/DP[$TI:0] >= $DP && INFO/TLOD >= $TLOD"
  # variant must be (nearly) absent from a normal that was sequenced deeply enough to see it
  if [ -n "$NI" ]; then EXPR+=" && FORMAT/DP[$NI:0] >= 10 && (FORMAT/AF[$NI:0] <= 0.03 || FORMAT/AF[$NI:0] = \".\")"; fi
  # POPAF is -log10(population AF): >= 3 means gnomAD AF <= 0.001
  [ "$SOMATIC_MODE" = tumor_only ] && EXPR+=" && INFO/POPAF >= 3"
  bcftools filter -i "$EXPR" -Oz -o "$OUT/$T.high_confidence.vcf.gz" "$OUT/$T.pass.vcf.gz"
  tabix -f -p vcf "$OUT/$T.high_confidence.vcf.gz"
  gatkj VariantsToTable -V "$OUT/$T.high_confidence.vcf.gz" -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER \
    -F TLOD -F NLOD -F ECNT -GF GT -GF AD -GF AF -GF DP -O "$OUT/$T.mutations.tsv"
  RAW=$(bcftools view -H "$OUT/$T.raw.vcf.gz" | wc -l); PASS=$(bcftools view -H "$OUT/$T.pass.vcf.gz" | wc -l)
  HC=$(bcftools view -H "$OUT/$T.high_confidence.vcf.gz" | wc -l)
  echo "$T [$SOMATIC_MODE]  raw=$RAW  PASS=$PASS  high_confidence=$HC  contamination: $(tail -1 "$C/$T.contamination.table" | cut -f2)"
  echo "  filter used: $EXPR"
done
