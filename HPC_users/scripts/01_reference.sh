#!/usr/bin/env bash
# Download reference genome, GATK resource bundle, somatic resources and gene models; build the BWA index.
# Run once per genome build (BWA index takes 2-3 hours). Needs internet - on HPC use a login/transfer node.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env core
mkdir -p "$RES_DIR" "$SOM_DIR" "$REF_DIR/annot"
dl() { [ -s "$2" ] && return 0; wget -c -q --show-progress -O "$2.part" "$1" && mv "$2.part" "$2" || { echo "WARN: download failed: $1"; rm -f "$2.part"; }; }

echo "== genome + germline resources ($BUILD)"
for u in "${GERMLINE_URLS[@]}"; do dl "$u" "$RES_DIR/$(basename "$u")"; done
echo "== somatic resources (PoN, gnomAD AF, common variants)"
for u in "${SOMATIC_URLS[@]}"; do dl "$u" "$SOM_DIR/$(basename "$u")"; done
for v in "$SOM_DIR"/*.vcf; do [ -e "$v" ] && bgzip -f "$v"; done       # hg19 files come uncompressed
for v in "$RES_DIR"/*.vcf.gz "$SOM_DIR"/*.vcf.gz; do [ -e "$v" ] && { [ -s "$v.tbi" ] || tabix -f -p vcf "$v"; }; done

echo "== gene models for CNVkit (UCSC refFlat)"
if [ ! -s "$REFFLAT" ]; then
  dl "https://hgdownload.soe.ucsc.edu/goldenPath/$BUILD/database/refFlat.txt.gz" "$REFFLAT.gz"
  [ -s "$REFFLAT.gz" ] && gunzip -f "$REFFLAT.gz"
fi

echo "== indexes"
[ -s "$REF_FASTA.fai" ] || samtools faidx "$REF_FASTA"
if [ ! -s "$REF_DICT" ]; then use_env gatk; gatk CreateSequenceDictionary -R "$REF_FASTA" -O "$REF_DICT"; use_env core; fi
[ -s "$REF_FASTA.bwt" ] || bwa index "$REF_FASTA"
echo "Reference ready: $REF_FASTA"
