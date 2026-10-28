#!/usr/bin/env bash
# Filter each sample's variants -> vcf/<ID>.filtered.vcf.gz (all, FILTER column) and vcf/<ID>.pass.vcf.gz (PASS only).
# FILTER_MODE=hard (default) uses separate SNP and indel hard filters. vqsr needs many samples (GATK: ~30+).
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
[ "$FILTER_MODE" = "vqsr" ] && echo "NOTE: VQSR on a single sample usually fails or under-performs. Prefer FILTER_MODE=hard, or run VQSR on the joint cohort (script 10)."
for S in $(samples_from_args "$@"); do
  filter_variants "$PROJ/vcf/$S.raw.vcf.gz" "$PROJ/vcf/$S"
done
