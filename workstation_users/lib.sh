#!/usr/bin/env bash
# Shared helpers + settings derived from config.sh. Sourced at the end of config.sh. Do not run directly.

# ---------------- helpers ----------------
use_env() {                                   # use_env core  -> activates conda env "<ENV_PREFIX>_core"
  local base="${CONDA_BASE:-}"
  [ -n "$base" ] || base="$(conda info --base 2>/dev/null || true)"
  if [ -z "$base" ]; then
    for d in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"; do
      if [ -d "$d" ]; then base="$d"; break; fi
    done
  fi
  [ -n "$base" ] || { echo "conda not found: set CONDA_BASE in config.sh"; return 1; }
  set +u
  # shellcheck disable=SC1091
  source "$base/etc/profile.d/conda.sh"
  conda activate "${ENV_PREFIX}_$1"
  set -u
}
read_list()          { grep -v '^[[:space:]]*#' "$1" | grep -v '^[[:space:]]*$' || true; }
samples_from_args()  { if [ "$#" -gt 0 ]; then printf '%s\n' "$@"; else read_list "$SAMPLE_LIST"; fi; }
tumors()             { read_list "$PAIRS_FILE" | awk '{print $1}'; }
tumors_from_args()   { if [ "$#" -gt 0 ]; then printf '%s\n' "$@"; else tumors; fi; }
normal_of()          { read_list "$PAIRS_FILE" | awk -v t="$1" '$1==t{print $2}'; }
has_annotator()      { printf '%s\n' "${ANNOTATORS[@]}" | grep -qx "$1"; }
gatkj()              { gatk --java-options "-Xmx${JAVA_MEM}g -Djava.io.tmpdir=${TMP_DIR}" "$@"; }

# ---------------- derived settings ----------------
export TMP_DIR="${TMP_DIR:-$PROJ/tmp}"; mkdir -p "$TMP_DIR"
JAVA_MEM=$(( MEM_GB > 8 ? MEM_GB - 4 : MEM_GB ))
HMM_THREADS=$(( CORES < 8 ? CORES : 8 ))
TARGETS_IL="$PROJ/targets/targets.interval_list"
GWAS_DIR="${GWAS_DIR:-$PROJ/gwas}"

if [ "$SEQ_TYPE" = "WES" ] && [ -z "$TARGETS" ]; then
  echo "SEQ_TYPE=WES needs TARGETS=<capture-kit BED> in config.sh"; return 1
fi
if [ -n "$TARGETS" ]; then INTERVAL_ARGS=(-L "$TARGETS" --interval-padding "$PAD"); else INTERVAL_ARGS=(); fi

BR="https://storage.googleapis.com/gcp-public-data--broad-references"
case "$BUILD" in
  hg38)
    RES_DIR="$REF_DIR/hg38"; SOM_DIR="$RES_DIR/somatic"
    G="$BR/hg38/v0"; S="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
    REF_FASTA="$RES_DIR/Homo_sapiens_assembly38.fasta"
    GERMLINE_URLS=("$G/Homo_sapiens_assembly38.fasta" "$G/Homo_sapiens_assembly38.fasta.fai" "$G/Homo_sapiens_assembly38.dict")
    for f in hapmap_3.3.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
             Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Homo_sapiens_assembly38.dbsnp138.vcf.gz; do
      GERMLINE_URLS+=("$G/$f" "$G/$f.tbi")
    done
    HAPMAP="$RES_DIR/hapmap_3.3.hg38.vcf.gz"; OMNI="$RES_DIR/1000G_omni2.5.hg38.vcf.gz"
    G1000="$RES_DIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    MILLS="$RES_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; DBSNP="$RES_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    SOMATIC_URLS=()
    for f in 1000g_pon.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz; do SOMATIC_URLS+=("$S/$f" "$S/$f.tbi"); done
    PON="$SOM_DIR/1000g_pon.hg38.vcf.gz"; GNOMAD="$SOM_DIR/af-only-gnomad.hg38.vcf.gz"; COMMON="$SOM_DIR/small_exac_common_3.hg38.vcf.gz"
    CHROMS=(); for i in {1..22} X Y; do CHROMS+=("chr$i"); done
    : "${SNPEFF_DB:=GRCh38.105}"; FUNCO_REF="hg38"; VEP_ASSEMBLY="GRCh38"; BSGENOME="BSgenome.Hsapiens.UCSC.hg38"; ANNOVAR_BUILD="hg38"
    ;;
  hg19)
    RES_DIR="$REF_DIR/hg19"; SOM_DIR="$RES_DIR/somatic"
    G="$BR/hg19/v0"; S="https://storage.googleapis.com/gatk-best-practices/somatic-b37"
    REF_FASTA="$RES_DIR/Homo_sapiens_assembly19.fasta"
    GERMLINE_URLS=("$G/Homo_sapiens_assembly19.fasta" "$G/Homo_sapiens_assembly19.fasta.fai" "$G/Homo_sapiens_assembly19.dict")
    # Check names in the bucket listing if a download fails (see README).
    for f in hapmap_3.3.b37.sites.vcf.gz 1000G_omni2.5.b37.sites.vcf.gz 1000G_phase1.snps.high_confidence.b37.vcf.gz \
             Mills_and_1000G_gold_standard.indels.b37.sites.vcf.gz dbsnp_138.b37.vcf.gz; do
      GERMLINE_URLS+=("$G/$f" "$G/$f.tbi")
    done
    HAPMAP="$RES_DIR/hapmap_3.3.b37.sites.vcf.gz"; OMNI="$RES_DIR/1000G_omni2.5.b37.sites.vcf.gz"
    G1000="$RES_DIR/1000G_phase1.snps.high_confidence.b37.vcf.gz"
    MILLS="$RES_DIR/Mills_and_1000G_gold_standard.indels.b37.sites.vcf.gz"; DBSNP="$RES_DIR/dbsnp_138.b37.vcf.gz"
    SOMATIC_URLS=("$S/Mutect2-WGS-panel-b37.vcf" "$S/af-only-gnomad.raw.sites.vcf" "$S/small_exac_common_3.vcf")
    PON="$SOM_DIR/Mutect2-WGS-panel-b37.vcf.gz"; GNOMAD="$SOM_DIR/af-only-gnomad.raw.sites.vcf.gz"; COMMON="$SOM_DIR/small_exac_common_3.vcf.gz"
    CHROMS=(); for i in {1..22} X Y; do CHROMS+=("$i"); done
    : "${SNPEFF_DB:=GRCh37.75}"; FUNCO_REF="b37"; VEP_ASSEMBLY="GRCh37"; BSGENOME="BSgenome.Hsapiens.UCSC.hg19"; ANNOVAR_BUILD="hg19"
    ;;
  *) echo "BUILD must be hg38 or hg19"; return 1 ;;
esac
REF_DICT="${REF_FASTA%.fasta}.dict"
REFFLAT="$REF_DIR/annot/refFlat.txt"

# Work units for scatter steps: whole target file for WES, one contig at a time for WGS
if [ -n "$TARGETS" ]; then UNITS=(all); else UNITS=("${CHROMS[@]}"); fi
set_unit() { if [ "$1" = all ]; then UNIT_L=(-L "$TARGETS"); else UNIT_L=(-L "$1"); fi; }

funco_dir() { ls -d "$FUNCO_DIR"/funcotator_dataSources.*"$1" 2>/dev/null | sort | tail -1; }   # $1 = g (germline) or s (somatic)

# ---------------- germline variant filtering ----------------
# filter_variants <raw.vcf.gz> <out_prefix>  ->  <prefix>.filtered.vcf.gz (all, FILTER column set) and <prefix>.pass.vcf.gz
filter_variants() {
  local in="$1" pre="$2"
  if [ "$FILTER_MODE" = "vqsr" ]; then
    gatkj VariantRecalibrator -R "$REF_FASTA" -V "$in" \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 "$OMNI" \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$G1000" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      -mode SNP -O "$pre.snps.recal" --tranches-file "$pre.snps.tranches"
    gatkj ApplyVQSR -R "$REF_FASTA" -V "$in" --recal-file "$pre.snps.recal" --tranches-file "$pre.snps.tranches" \
      -mode SNP --truth-sensitivity-filter-level "$VQSR_SNP_SENS" -O "$pre.snps_recal.vcf.gz"
    gatkj VariantRecalibrator -R "$REF_FASTA" -V "$pre.snps_recal.vcf.gz" \
      --resource:mills,known=false,training=true,truth=true,prior=12.0 "$MILLS" \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
      -an QD -an ReadPosRankSum -an FS -an SOR -mode INDEL --max-gaussians 4 \
      -O "$pre.indels.recal" --tranches-file "$pre.indels.tranches"
    gatkj ApplyVQSR -R "$REF_FASTA" -V "$pre.snps_recal.vcf.gz" --recal-file "$pre.indels.recal" \
      --tranches-file "$pre.indels.tranches" -mode INDEL --truth-sensitivity-filter-level "$VQSR_INDEL_SENS" \
      -O "$pre.filtered.vcf.gz"
  else
    local dp=()
    if [ "$SEQ_TYPE" = "WES" ]; then
      dp=(--filter-expression "DP < $WES_MIN_DP" --filter-name LowDepth --filter-expression "DP > $WES_MAX_DP" --filter-name HighDepth)
    fi
    gatkj SelectVariants -R "$REF_FASTA" -V "$in" --select-type-to-include SNP -O "$pre.snp.vcf.gz"
    gatkj SelectVariants -R "$REF_FASTA" -V "$in" --select-type-to-include INDEL --select-type-to-include MIXED -O "$pre.indel.vcf.gz"
    gatkj VariantFiltration -R "$REF_FASTA" -V "$pre.snp.vcf.gz" \
      --filter-expression "QD < 2.0" --filter-name LowQD \
      --filter-expression "FS > 60.0" --filter-name HighFS \
      --filter-expression "MQ < 40.0" --filter-name LowMQ \
      --filter-expression "SOR > 3.0" --filter-name HighSOR \
      --filter-expression "MQRankSum < -12.5" --filter-name LowMQRankSum \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name LowReadPosRankSum \
      "${dp[@]}" -O "$pre.snp.f.vcf.gz"
    gatkj VariantFiltration -R "$REF_FASTA" -V "$pre.indel.vcf.gz" \
      --filter-expression "QD < 2.0" --filter-name LowQD \
      --filter-expression "FS > 200.0" --filter-name HighFS \
      --filter-expression "SOR > 10.0" --filter-name HighSOR \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name LowReadPosRankSum \
      "${dp[@]}" -O "$pre.indel.f.vcf.gz"
    gatkj MergeVcfs -I "$pre.snp.f.vcf.gz" -I "$pre.indel.f.vcf.gz" -O "$pre.filtered.vcf.gz"
    rm -f "$pre".snp.vcf.gz* "$pre".indel.vcf.gz* "$pre".snp.f.vcf.gz* "$pre".indel.f.vcf.gz*
  fi
  bcftools view -f PASS -Oz -o "$pre.pass.vcf.gz" "$pre.filtered.vcf.gz"
  tabix -f -p vcf "$pre.pass.vcf.gz"
}

# ---------------- frequency tags (used by disease/patient-specific step) ----------------
# group_freq_annot <in.vcf.gz> <samples_file> <TagPrefix> <out.vcf.gz>  -> INFO/<Prefix>_AF, _AC, _AN
group_freq_annot() {
  local in="$1" grp="$2" pre="$3" out="$4" t; t="$(mktemp -d "$TMP_DIR/gf.XXXXXX")"
  bcftools view -S "$grp" --force-samples -Ou "$in" | bcftools +fill-tags -Ou -- -t AF,AC,AN \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\n' | bgzip -c > "$t/f.tsv.gz"
  tabix -s1 -b2 -e2 "$t/f.tsv.gz"
  {
    echo "##INFO=<ID=${pre}_AF,Number=A,Type=Float,Description=\"Allele frequency in ${pre} samples\">"
    echo "##INFO=<ID=${pre}_AC,Number=A,Type=Integer,Description=\"Allele count in ${pre} samples\">"
    echo "##INFO=<ID=${pre}_AN,Number=1,Type=Integer,Description=\"Allele number in ${pre} samples\">"
  } > "$t/h.txt"
  bcftools annotate -a "$t/f.tsv.gz" -h "$t/h.txt" -c "CHROM,POS,REF,ALT,INFO/${pre}_AF,INFO/${pre}_AC,INFO/${pre}_AN" -Oz -o "$out" "$in"
  tabix -f -p vcf "$out"; rm -rf "$t"
}
# add_gnomad_af <in.vcf.gz> <out.vcf.gz>  -> INFO/gnomAD_AF (looked up only at the cohort's sites)
add_gnomad_af() {
  local in="$1" out="$2" t; t="$(mktemp -d "$TMP_DIR/gn.XXXXXX")"
  bcftools query -f '%CHROM\t%POS\n' "$in" | uniq | bgzip -c > "$t/sites.tsv.gz"
  tabix -s1 -b2 -e2 "$t/sites.tsv.gz"
  bcftools view -R "$t/sites.tsv.gz" -Ou "$GNOMAD" | bcftools norm -m -any -Ou \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | bgzip -c > "$t/g.tsv.gz"
  tabix -s1 -b2 -e2 "$t/g.tsv.gz"
  echo '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD allele frequency">' > "$t/h.txt"
  bcftools annotate -a "$t/g.tsv.gz" -h "$t/h.txt" -c CHROM,POS,REF,ALT,INFO/gnomAD_AF -Oz -o "$out" "$in"
  tabix -f -p vcf "$out"; rm -rf "$t"
}
