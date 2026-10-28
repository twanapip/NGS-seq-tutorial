#!/usr/bin/env bash
# Download databases for the annotators listed in ANNOTATORS (config.sh). Big downloads - run once.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
if has_annotator snpeff; then
  use_env snpeff; snpEff download -v "$SNPEFF_DB"
fi
if has_annotator funcotator; then
  use_env gatk; mkdir -p "$FUNCO_DIR"; cd "$FUNCO_DIR"
  gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download
  gatk FuncotatorDataSourceDownloader --somatic  --validate-integrity --extract-after-download
fi
if has_annotator vep; then
  use_env vep; mkdir -p "$VEP_CACHE/Plugins"
  vep_install -a c -s homo_sapiens -y "$VEP_ASSEMBLY" --CACHE_VERSION "$VEP_CACHE_VERSION" --CACHEDIR "$VEP_CACHE" --NO_UPDATE
  for p in dbNSFP CADD REVEL; do   # plugin code only; their data files are separate downloads (see README)
    wget -q -nc -P "$VEP_CACHE/Plugins" "https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/$VEP_CACHE_VERSION/$p.pm" || true
  done
fi
if has_annotator annovar; then
  if [ -z "$ANNOVAR_DIR" ] || [ ! -f "$ANNOVAR_DIR/annotate_variation.pl" ]; then
    echo "ANNOVAR needs manual registration: https://www.openbioinformatics.org/annovar/annovar_download_form.php"
    echo "Extract it, then set ANNOVAR_DIR in config.sh to the folder holding annotate_variation.pl"
  else
    for db in refGene gnomad312_genome 1000g2015aug clinvar_20220320 dbnsfp42c; do
      perl "$ANNOVAR_DIR/annotate_variation.pl" -buildver "$ANNOVAR_BUILD" -downdb -webfrom annovar "$db" "$ANNOVAR_DIR/humandb/"
    done
  fi
fi
