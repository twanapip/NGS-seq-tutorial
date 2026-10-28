#!/usr/bin/env bash
# Copy number variants with CNVkit.  Usage: bash scripts/17_cnv.sh [germline|tumor]   (default: CNV_MODE in config.sh)
#  germline: tests = samples.txt ; reference = normals.txt (or all samples pooled if normals.txt is empty)
#  tumor   : tests = tumors in pairs.tsv ; reference = normals.txt (flat reference if empty)
# Output: cnv/<mode>/results/*  (.cnr bins, .cns segments, .call.cns calls) plus filtered tables, VCF/BED/SEG, plots.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
MODE="${1:-$CNV_MODE}"; use_env cnv
OUT="$PROJ/cnv/$MODE"; R="$OUT/results"; mkdir -p "$R"
[ -s "$REFFLAT" ] || { echo "Missing $REFFLAT - run scripts/01_reference.sh"; exit 1; }
bam() { echo "$PROJ/bam/$1.recal.bam"; }
mapfile -t NORMALS < <(read_list "$NORMALS_FILE")
case "$MODE" in
  germline) mapfile -t TESTS < <(read_list "$SAMPLE_LIST"); THR="$CNV_THRESH_GERMLINE"
            [ "${#NORMALS[@]}" -gt 0 ] || NORMALS=("${TESTS[@]}") ;;
  tumor)    mapfile -t TESTS < <(tumors | sort -u); THR="$CNV_THRESH_TUMOR" ;;
  *) echo "mode must be germline or tumor"; exit 1 ;;
esac
TB=(); for t in "${TESTS[@]}"; do TB+=("$(bam "$t")"); done
NB=(); for n in "${NORMALS[@]}"; do NB+=("$(bam "$n")"); done

ARGS=(batch "${TB[@]}" --fasta "$REF_FASTA" --annotate "$REFFLAT" --output-reference "$OUT/reference.cnn"
      --output-dir "$R" --segment-method cbs --drop-low-coverage --scatter --diagram -p "$CORES")
if [ "$SEQ_TYPE" = "WES" ]; then
  ARGS+=(--method hybrid --targets "$TARGETS")
else
  cnvkit.py access "$REF_FASTA" -s 10000 -o "$OUT/access.bed"
  ARGS+=(--method wgs --access "$OUT/access.bed")
fi
ARGS+=(--normal "${NB[@]}")            # no BAMs after --normal would build a flat reference
cnvkit.py "${ARGS[@]}"

PUR=(); [ -n "$PURITY" ] && PUR=(--purity "$PURITY")
for t in "${TESTS[@]}"; do
  P="$R/$t.recal"
  cnvkit.py call "$P.cns" -m threshold "-t=$THR" "${PUR[@]}" -o "$P.thr.call.cns"
  # keep real, supported, large-enough events (columns found by name, not position)
  awk -F'\t' -v min="$CNV_MIN_SIZE" -v pr="$CNV_MIN_PROBES" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;print;next} $h["cn"]!=2 && ($3-$2)>=min && $h["probes"]>=pr' \
    "$P.thr.call.cns" > "$OUT/$t.filtered.cns"
  cnvkit.py export vcf "$P.thr.call.cns" -o "$OUT/$t.cnv.vcf"
  cnvkit.py export bed "$P.thr.call.cns" -o "$OUT/$t.cnv.bed"
  cnvkit.py export seg "$P.thr.call.cns" -o "$OUT/$t.cnv.seg"
  cnvkit.py metrics "$P.cnr" -s "$P.cns" -o "$OUT/$t.metrics.txt"
  cnvkit.py genemetrics "$P.cnr" -s "$P.cns" -t 0.2 -m 3 -o "$OUT/$t.genemetrics.txt"
  for c in "${CHROMS[@]}"; do
    cnvkit.py scatter "$P.cnr" -s "$P.cns" -c "$c" -o "$OUT/$t.$c.scatter.pdf" --title "$t $c" || true
  done
done
cnvkit.py sex "$R"/*.cnr -o "$OUT/sex_check.txt" || true          # reported sex vs expected is a quick sample-swap check
[ "${#TESTS[@]}" -gt 1 ] && cnvkit.py heatmap "$R"/*.recal.cns -o "$OUT/heatmap.pdf" || true

if [ "$MODE" = germline ] && [ "$BUILD" = hg38 ]; then          # drop CNVs that are common in healthy people (DGV)
  DGV="$REF_DIR/annot/dgv_common_cnvs.bed"
  if [ ! -s "$DGV" ]; then
    wget -q -O "$REF_DIR/annot/dgv.txt" "$DGV_URL" && tail -n +2 "$REF_DIR/annot/dgv.txt" \
      | awk -F'\t' '$2!=""{print "chr"$2"\t"$3"\t"$4"\t"$1}' | sort -k1,1 -k2,2n > "$DGV" || echo "WARN: DGV download failed - skipping common-CNV filter"
  fi
  if [ -s "$DGV" ]; then
    for t in "${TESTS[@]}"; do
      tail -n +2 "$OUT/$t.filtered.cns" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
        | bedtools intersect -a - -b "$DGV" -v -f 0.5 > "$OUT/$t.rare.bed"
    done
  fi
fi

if [ "$MODE" = tumor ] && [ "$CNV_CLINVAR" = 1 ]; then              # overlap with ClinVar CNV-type records
  CV="$REF_DIR/annot/clinvar.vcf.gz"; DIR=$([ "$BUILD" = hg38 ] && echo vcf_GRCh38 || echo vcf_GRCh37)
  [ -s "$CV" ] || wget -q -O "$CV" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/$DIR/clinvar.vcf.gz"
  tabix -f -p vcf "$CV"
  PFX=$([ "$BUILD" = hg38 ] && echo chr || echo "")
  bcftools view -i 'INFO/CLNVC="Deletion" || INFO/CLNVC="Duplication" || INFO/CLNVC="copy_number_gain" || INFO/CLNVC="copy_number_loss"' "$CV" \
    | bcftools query -f '%CHROM\t%POS0\t%END\t%ID\t%INFO/CLNVC\t%INFO/CLNSIG\t%INFO/GENEINFO\n' | awk -v p="$PFX" 'BEGIN{OFS="\t"}{$1=p $1; print}' > "$OUT/clinvar_cnvs.bed"
  for t in "${TESTS[@]}"; do
    tail -n +2 "$OUT/$t.filtered.cns" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' \
      | bedtools intersect -a - -b "$OUT/clinvar_cnvs.bed" -wa -wb -f 0.5 > "$OUT/$t.clinvar_overlaps.bed"
  done
fi
echo "CNV results: $OUT"
