#!/usr/bin/env bash
# Build PLINK files (bed/bim/fam) for GWAS.
#   GWAS_VCF set                    -> your own joint-called cohort VCF
#   GWAS_VCF empty, GWAS_DATA=1000g -> downloads 1000 Genomes phase 3 (GRCh37, 2,504 people, ~15 GB) for practice
# Needs hundreds of samples for a real GWAS; tiny cohorts are for learning the workflow only.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gwas
mkdir -p "$GWAS_DIR"/{data,results,plots}; cd "$GWAS_DIR/data"

clean_to_plink() {   # clean_to_plink <in.vcf.gz> <out_prefix> <threads>: split multi-allelics, drop duplicates, keep bi-allelic, convert
  local in="$1" out="$2" th="$3"
  bcftools norm -m -any --threads "$th" "$in" | bcftools norm -d all --threads "$th" \
    | bcftools view -m2 -M2 -v snps,indels --threads "$th" -Oz -o "$out.clean.vcf.gz"
  tabix -f -p vcf "$out.clean.vcf.gz"
  plink --vcf "$out.clean.vcf.gz" --threads "$th" --double-id --biallelic-only strict \
        --set-missing-var-ids '@:#_$1_$2' --make-bed --out "$out" --silent
  rm -f "$out.clean.vcf.gz" "$out.clean.vcf.gz.tbi"
}

if [ -n "$GWAS_VCF" ]; then
  clean_to_plink "$GWAS_VCF" "$GWAS_DIR/data/cohort_raw" "$CORES"
elif [ "$GWAS_DATA" = "1000g" ]; then
  BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
  vcf_name() { if [ "$1" = X ]; then echo "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"; else echo "ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"; fi; }
  for C in {1..22} X; do
    f="$(vcf_name "$C")"
    [ -s "$f" ] || { wget -q -c "$BASE/$f"; wget -q -c "$BASE/$f.tbi"; }
  done
  wget -q -nc "$BASE/integrated_call_samples_v3.20130502.ALL.panel"
  JOBS=$(( CORES / 4 )); [ "$JOBS" -ge 1 ] || JOBS=1
  export -f clean_to_plink vcf_name
  printf '%s\n' {1..22} X | xargs -P "$JOBS" -I{} bash -c 'f=$(vcf_name {}); [ -s "1kg_chr{}.bed" ] || clean_to_plink "$f" "1kg_chr{}" 2'
  : > merge_list.txt; for C in {2..22} X; do echo "1kg_chr$C" >> merge_list.txt; done
  plink --bfile 1kg_chr1 --merge-list merge_list.txt --make-bed --out "$GWAS_DIR/data/cohort_raw" --silent
else
  echo "Set GWAS_VCF or GWAS_DATA=1000g in config.sh"; exit 1
fi
echo "Variants: $(wc -l < "$GWAS_DIR/data/cohort_raw.bim")   Samples: $(wc -l < "$GWAS_DIR/data/cohort_raw.fam")"
