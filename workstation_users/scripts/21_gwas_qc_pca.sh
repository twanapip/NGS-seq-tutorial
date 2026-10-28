#!/usr/bin/env bash
# GWAS quality control (samples, sex, variants) and PCA for ancestry. Output: gwas/results/clean.* , pca.eigenvec, plots/
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gwas
cd "$GWAS_DIR/results"
RAW="$GWAS_DIR/data/cohort_raw"
P="plink --silent"
$P --bfile "$RAW" --missing --out qc_missing
$P --bfile "$RAW" --mind "$GWAS_MIND" --make-bed --out s1
CUR=s1
if awk '$1==23||$1=="X"{f=1} END{exit !f}' s1.bim; then          # sex check needs X-chromosome data
  $P --bfile s1 --check-sex --out sex_check
  $P --bfile s1 --impute-sex --make-bed --out s2; CUR=s2
fi
if [ "$HET_FILTER" = 1 ]; then
  $P --bfile "$CUR" --het --out het
  awk 'NR==FNR{if(FNR>1){s+=$6;ss+=$6*$6;n++};next} FNR>1{m=s/n;sd=sqrt(ss/n-m*m); if(($6-m)>3*sd||(m-$6)>3*sd)print $1,$2}' het.het het.het > het_outliers.txt
  $P --bfile "$CUR" --remove het_outliers.txt --make-bed --out s3; CUR=s3
fi
if [ "$RELATEDNESS_FILTER" = 1 ]; then
  $P --bfile "$CUR" --indep-pairwise 50 5 0.2 --out rel_prune
  $P --bfile "$CUR" --extract rel_prune.prune.in --rel-cutoff 0.125 --out rel
  $P --bfile "$CUR" --keep rel.rel.id --make-bed --out s4; CUR=s4
fi
$P --bfile "$CUR" --geno "$GWAS_GENO" --make-bed --out v1
$P --bfile v1 --maf "$GWAS_MAF" --make-bed --out v2
$P --bfile v2 --hwe "$GWAS_HWE" --make-bed --out clean
$P --bfile clean --missing --freq --hardy --out final_qc_stats
{
  echo "Samples: raw $(wc -l < "$RAW.fam") -> final $(wc -l < clean.fam)"
  echo "Variants: raw $(wc -l < "$RAW.bim") -> final $(wc -l < clean.bim)"
} | tee "$GWAS_DIR/qc_summary.txt"

# PCA on LD-pruned variants (so ancestry, not local LD blocks, drives the components)
$P --bfile clean --indep-pairwise 50 5 0.2 --out pca_pruning
$P --bfile clean --extract pca_pruning.prune.in --pca 20 --out pca
PANEL="$POP_PANEL"; [ -n "$PANEL" ] || PANEL="$GWAS_DIR/data/integrated_call_samples_v3.20130502.ALL.panel"
export GWAS_DIR PANEL
cat > plot_pca.R <<'EOF_R'
library(ggplot2)
setwd(file.path(Sys.getenv("GWAS_DIR"), "results"))
pca <- read.table("pca.eigenvec", header = FALSE); colnames(pca) <- c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2)))
ev <- read.table("pca.eigenval")$V1; var <- 100 * ev / sum(ev)
panel <- Sys.getenv("PANEL")
if (nzchar(panel) && file.exists(panel)) {
  info <- read.table(panel, header = TRUE, stringsAsFactors = FALSE)
  pca <- merge(pca, info[, c("sample", "super_pop")], by.x = "IID", by.y = "sample", all.x = TRUE)
} else pca$super_pop <- "all samples"
for (pair in list(c(1, 2), c(3, 4))) {
  a <- paste0("PC", pair[1]); b <- paste0("PC", pair[2])
  g <- ggplot(pca, aes(.data[[a]], .data[[b]], colour = super_pop)) + geom_point(alpha = 0.6) + theme_minimal() +
    labs(x = sprintf("%s (%.2f%%)", a, var[pair[1]]), y = sprintf("%s (%.2f%%)", b, var[pair[2]]), colour = "Group", title = "Population structure (PCA)")
  ggsave(file.path("..", "plots", sprintf("pca_%s_%s.png", a, b)), g, width = 8, height = 6, dpi = 120)
}
sc <- data.frame(PC = seq_along(var), Variance = var)
ggsave(file.path("..", "plots", "pca_scree.png"),
  ggplot(sc, aes(PC, Variance)) + geom_col(fill = "steelblue") + geom_line() + geom_point() + theme_minimal() + labs(y = "% variance", title = "Scree plot"),
  width = 7, height = 5, dpi = 120)
EOF_R
Rscript plot_pca.R
echo "QC + PCA done. Look at plots/pca_scree.png to choose N_PCS (elbow; usually 5-10)."
