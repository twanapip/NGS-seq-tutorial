#!/usr/bin/env bash
# Variant statistics, one-page QC summary, tables, BED and coverage tracks for IGV.  Output: qc/<ID>.qc_summary.txt
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gatk
M="$PROJ/qc/metrics"
for S in $(samples_from_args "$@"); do
  V="$PROJ/vcf/$S.pass.vcf.gz"
  bcftools stats "$V" > "$M/$S.variant_stats.txt"
  NSNP=$(bcftools view -H -v snps "$V" | wc -l); NIND=$(bcftools view -H -v indels "$V" | wc -l)
  TITV=$(grep '^TSTV' "$M/$S.variant_stats.txt" | cut -f5 | head -1)
  HETHOM=$(bcftools query -f '[%GT\n]' "$V" | awk '$0~/^(0[\/|]1|1[\/|]0|1[\/|]2|2[\/|]1)$/{h++} $0~/^(1[\/|]1|2[\/|]2)$/{o++} END{if(o>0)printf "%.2f",h/o; else print "NA"}')
  TOT=$(bcftools view -H "$V" | wc -l); KNOWN=$(bcftools view -H "$V" | awk '$3!="."' | wc -l)
  DBPCT=$(awk -v k="$KNOWN" -v t="$TOT" 'BEGIN{if(t>0)printf "%.1f",100*k/t; else print "NA"}')
  MAP=$(awk -F'\t' '/^CATEGORY/{for(i=1;i<=NF;i++)h[$i]=i;next} $1=="PAIR"{printf "%.2f",100*$h["PCT_PF_READS_ALIGNED"]}' "$M/$S.alignment_summary.txt")
  DUP=$(awk -F'\t' '/^LIBRARY/{for(i=1;i<=NF;i++)h[$i]=i;next} h["PERCENT_DUPLICATION"]&&NF>5{printf "%.2f",100*$h["PERCENT_DUPLICATION"];exit}' "$M/$S.dup_metrics.txt")
  INS=$(awk '/^MEDIAN_INSERT_SIZE/{getline;print $1;exit}' "$M/$S.insert_size.txt")
  if [ "$SEQ_TYPE" = "WES" ]; then
    COV=$(awk -F'\t' '/^BAIT_SET/{for(i=1;i<=NF;i++)h[$i]=i;next} h["MEAN_TARGET_COVERAGE"]&&NF>5{print $h["MEAN_TARGET_COVERAGE"];exit}' "$M/$S.hs_metrics.txt")
    ONT=$(awk -F'\t' '/^BAIT_SET/{for(i=1;i<=NF;i++)h[$i]=i;next} h["PCT_SELECTED_BASES"]&&NF>5{printf "%.1f",100*$h["PCT_SELECTED_BASES"];exit}' "$M/$S.hs_metrics.txt")
    P20=$(awk -F'\t' '/^BAIT_SET/{for(i=1;i<=NF;i++)h[$i]=i;next} h["PCT_TARGET_BASES_20X"]&&NF>5{printf "%.1f",100*$h["PCT_TARGET_BASES_20X"];exit}' "$M/$S.hs_metrics.txt")
  else
    COV=$(awk -F'\t' '/^GENOME_TERRITORY/{for(i=1;i<=NF;i++)h[$i]=i;next} h["MEAN_COVERAGE"]&&NF>5{print $h["MEAN_COVERAGE"];exit}' "$M/$S.wgs_metrics.txt")
  fi
  {
    echo "QC summary: $S ($SEQ_TYPE, $BUILD)"
    echo "metric | value | benchmark"
    echo "Reads aligned % | $MAP | >95"
    echo "Duplicate rate % | $DUP | <30"
    echo "Median insert size bp | $INS | ~300-500 (library dependent)"
    if [ "$SEQ_TYPE" = "WES" ]; then
      echo "Mean target coverage | $COV | >100 clinical, >50 research"
      echo "Reads on target % | $ONT | >80"
      echo "Targets at >=20x % | $P20 | >95"
    else
      echo "Mean coverage | $COV | >=30"
    fi
    echo "SNPs (PASS) | $NSNP | WGS 4-5 million; exome 20-25 thousand"
    echo "Indels (PASS) | $NIND | WGS 0.5-0.8 million"
    echo "Ti/Tv | $TITV | WGS ~2.0-2.1; exome ~2.8-3.3"
    echo "Het/Hom | $HETHOM | 1.5-2.0 (population dependent)"
    echo "Variants in dbSNP % | $DBPCT | >95 for common populations"
  } > "$PROJ/qc/$S.qc_summary.txt"
  cat "$PROJ/qc/$S.qc_summary.txt"
  gatkj VariantsToTable -V "$V" -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F TYPE -GF GT -GF DP -GF GQ -O "$PROJ/vcf/$S.variants_table.tsv"
  bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$V" > "$PROJ/vcf/$S.variants.bed"
  if [ "$MAKE_BEDGRAPH" = 1 ]; then
    bedtools genomecov -ibam "$PROJ/bam/$S.recal.bam" -bg > "$PROJ/bam/$S.coverage.bedgraph"
  fi
done
