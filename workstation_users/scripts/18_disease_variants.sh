#!/usr/bin/env bash
# Patient-specific and disease-specific (case vs control) candidate variants from the annotated cohort VCF.
# Needs: cohort VCF annotated with SnpEff (ANN field)  -> run scripts/10 then scripts/12 on cohort/cohort.pass.vcf.gz first.
# Usage: bash scripts/18_disease_variants.sh [annotated.vcf.gz]
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
IN="${1:-$PROJ/annotation/cohort.pass/snpeff.vcf.gz}"
[ -s "$IN" ] || { echo "Annotated cohort VCF not found: $IN"; exit 1; }
D="$PROJ/cohort/disease"; mkdir -p "$D/patient_specific" "$D/disease_specific"
for f in "$CASES_FILE" "$CONTROLS_FILE"; do [ -n "$(read_list "$f")" ] || { echo "Fill in $f"; exit 1; }; done
read_list "$CASES_FILE" > "$D/cases.txt"; read_list "$CONTROLS_FILE" > "$D/controls.txt"

# 1) one ALT allele per line, then add gnomAD, control and case frequencies
bcftools norm -m -any -Oz -o "$D/biallelic.vcf.gz" "$IN"; tabix -f -p vcf "$D/biallelic.vcf.gz"
add_gnomad_af "$D/biallelic.vcf.gz" "$D/step1.vcf.gz"
group_freq_annot "$D/step1.vcf.gz" "$D/controls.txt" Control "$D/step2.vcf.gz"
group_freq_annot "$D/step2.vcf.gz" "$D/cases.txt" Case "$D/annotated_freq.vcf.gz"
rm -f "$D"/step1.vcf.gz* "$D"/step2.vcf.gz*

CTRL="(INFO/Control_AF < $DISEASE_CONTROL_AF_MAX || INFO/Control_AF = \".\")"
GNOM="(INFO/gnomAD_AF < $DISEASE_GNOMAD_AF_MAX || INFO/gnomAD_AF = \".\")"
IMP='(INFO/ANN ~ "HIGH" || INFO/ANN ~ "MODERATE")'
genes() { awk -F, '{split($1,a,"|"); print a[4]}'; }

# 2) patient-specific: carried by the patient, rare in controls and gnomAD, HIGH/MODERATE impact
for P in $(cat "$D/cases.txt"); do
  O="$D/patient_specific/$P"; mkdir -p "$O"
  bcftools view -s "$P" -i 'GT="alt"' "$D/annotated_freq.vcf.gz" \
    | bcftools view -i "$CTRL && $GNOM && $IMP" -Oz -o "$O/candidates.vcf.gz"
  tabix -f -p vcf "$O/candidates.vcf.gz"
  { echo -e "CHROM\tPOS\tREF\tALT\tGENOTYPE\tgnomAD_AF\tControl_AF\tANN"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%INFO/gnomAD_AF\t%INFO/Control_AF\t%INFO/ANN\n' "$O/candidates.vcf.gz"; } > "$O/candidates_table.tsv"
  # genes hit by 2+ different heterozygous variants = compound-heterozygous candidates (phase is unknown: confirm with parents/reads)
  bcftools view -i 'GT="het"' "$O/candidates.vcf.gz" | bcftools query -f '%INFO/ANN\n' | genes | sort | uniq -c | awk '$1>=2' > "$O/compound_het_genes.txt"
done

# 3) disease-specific: more frequent in cases, seen in >= N cases, rare elsewhere, damaging
S="$D/disease_specific"
bcftools view -i "INFO/Case_AF > INFO/Control_AF && INFO/Case_AC >= $DISEASE_MIN_CASE_AC && $CTRL && $GNOM && $IMP" \
  -Oz -o "$S/disease_enriched.vcf.gz" "$D/annotated_freq.vcf.gz"
tabix -f -p vcf "$S/disease_enriched.vcf.gz"
{ echo -e "CHROM\tPOS\tREF\tALT\tCase_AC\tCase_AF\tControl_AC\tControl_AF\tgnomAD_AF\tANN"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/Case_AC\t%INFO/Case_AF\t%INFO/Control_AC\t%INFO/Control_AF\t%INFO/gnomAD_AF\t%INFO/ANN\n' \
    "$S/disease_enriched.vcf.gz"; } > "$S/disease_enriched_table.tsv"
bcftools query -f '%INFO/ANN\n' "$S/disease_enriched.vcf.gz" | genes | sort | uniq -c | sort -rn > "$S/recurrent_genes.txt"
# variant-type subsets
bcftools view -i 'INFO/ANN ~ "HIGH"' -Oz -o "$S/loss_of_function.vcf.gz" "$S/disease_enriched.vcf.gz"
bcftools view -i 'INFO/ANN ~ "missense_variant"' -Oz -o "$S/missense.vcf.gz" "$S/disease_enriched.vcf.gz"
bcftools view -i 'INFO/ANN ~ "splice"' -Oz -o "$S/splice.vcf.gz" "$S/disease_enriched.vcf.gz"
echo "Disease-enriched variants: $(bcftools view -H "$S/disease_enriched.vcf.gz" | wc -l)   (top genes: $S/recurrent_genes.txt)"
