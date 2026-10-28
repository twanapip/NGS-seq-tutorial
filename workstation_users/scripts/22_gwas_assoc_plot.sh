#!/usr/bin/env bash
# Association test (linear or logistic, adjusted for sex + PCs), significant hits, Manhattan + QQ plots, lambda.
# Needs PHENO_FILE (no header: FID IID phenotype) or SIMULATE_PHENO=1 for a practice trait.
set -euo pipefail
source "$(dirname "$0")/../config.sh"
use_env gwas
cd "$GWAS_DIR/results"
export GWAS_DIR N_PCS COVAR_EXTRA PHENO_FILE SIMULATE_PHENO
cat > make_covariates.R <<'EOF_R'
setwd(file.path(Sys.getenv("GWAS_DIR"), "results")); npc <- as.integer(Sys.getenv("N_PCS"))
fam <- read.table("clean.fam"); colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
pca <- read.table("pca.eigenvec"); colnames(pca) <- c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2)))
cov <- merge(fam[, c("FID", "IID", "SEX")], pca[, 1:(2 + npc)], by = c("FID", "IID"))
extra <- Sys.getenv("COVAR_EXTRA")
if (nzchar(extra)) cov <- merge(cov, read.table(extra, header = TRUE), by = c("FID", "IID"))
write.table(cov, "covariates.txt", sep = "\t", quote = FALSE, row.names = FALSE)
if (Sys.getenv("SIMULATE_PHENO") == "1" && !nzchar(Sys.getenv("PHENO_FILE"))) {   # practice trait: ancestry + sex + noise, NO real genetics
  set.seed(12345)
  y <- 10 + 0.3 * cov$PC1 + 0.2 * cov$PC2 + ifelse(cov$SEX == 1, 0, 0.5) + rnorm(nrow(cov))
  write.table(data.frame(cov$FID, cov$IID, scale(y)[, 1]), "simulated_phenotype.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
EOF_R
Rscript make_covariates.R
PH="$PHENO_FILE"; [ -n "$PH" ] || { [ -s simulated_phenotype.txt ] && PH="$PWD/simulated_phenotype.txt"; }
[ -n "$PH" ] || { echo "Provide PHENO_FILE (or SIMULATE_PHENO=1 for practice) in config.sh"; exit 1; }
COVN=$(head -1 covariates.txt | cut -f3- | tr '\t' ',')
if [ "$TRAIT_TYPE" = binary ]; then MODEL="--logistic hide-covar"; RES=gwas.assoc.logistic; else MODEL="--linear hide-covar"; RES=gwas.assoc.linear; fi
# shellcheck disable=SC2086
plink --silent --bfile clean --pheno "$PH" --covar covariates.txt --covar-name "$COVN" $MODEL --out gwas
awk 'NR==1 || ($5=="ADD" && $9!="NA" && $9<5e-8)' "$RES" > significant_snps.txt
awk 'NR==1 || ($5=="ADD" && $9!="NA" && $9<1e-5)' "$RES" > suggestive_snps.txt
(head -1 "$RES"; awk 'NR>1 && $5=="ADD" && $9!="NA"' "$RES" | sort -k9,9g | head -20) > top_associations.txt
echo "Genome-wide significant (p<5e-8): $(($(wc -l < significant_snps.txt) - 1))   suggestive (p<1e-5): $(($(wc -l < suggestive_snps.txt) - 1))"
export RES
cat > plot_gwas.R <<'EOF_R'
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
setwd(file.path(Sys.getenv("GWAS_DIR"), "results"))
g <- read.table(Sys.getenv("RES"), header = TRUE, stringsAsFactors = FALSE) %>%
  filter(TEST == "ADD", !is.na(P), P > 0, P <= 1) %>% mutate(lp = -log10(P))
chisq <- qchisq(1 - g$P, 1); lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
# thin the unremarkable points so the plot stays light (all p < 0.01 are kept)
set.seed(1); gp <- g %>% filter(P < 0.01 | runif(n()) < 0.05) %>% arrange(CHR, BP)
off <- gp %>% group_by(CHR) %>% summarise(len = max(BP)) %>% mutate(start = cumsum(lag(len, default = 0)))
gp <- left_join(gp, off, by = "CHR") %>% mutate(pos = BP + start)
ax <- gp %>% group_by(CHR) %>% summarise(mid = mean(pos))
ggsave("../plots/manhattan_plot.png", width = 16, height = 7, dpi = 110,
  ggplot(gp, aes(pos, lp, colour = factor(CHR %% 2))) + geom_point(size = 0.6, alpha = 0.7) +
    geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 5, colour = "blue", linetype = "dashed") +
    scale_colour_manual(values = c("steelblue", "darkblue"), guide = "none") +
    scale_x_continuous(breaks = ax$mid, labels = ax$CHR) + theme_minimal() +
    labs(x = "Chromosome", y = expression(-log[10](P)), title = "Manhattan plot", caption = "red: p = 5e-8 (genome-wide)   blue: p = 1e-5 (suggestive)"))
q <- g %>% arrange(P) %>% mutate(expected = -log10(ppoints(n())), observed = -log10(P))
q <- q %>% filter(row_number() <= 20000 | row_number() %% 20 == 0)
ggsave("../plots/qq_plot.png", width = 7, height = 7, dpi = 110,
  ggplot(q, aes(expected, observed)) + geom_point(alpha = 0.5, colour = "steelblue") + geom_abline(colour = "red", linetype = "dashed") +
    theme_minimal() + labs(title = sprintf("QQ plot (lambda = %.3f)", lambda), x = "Expected -log10(P)", y = "Observed -log10(P)"))
cat(sprintf("Variants tested: %d\nGenomic inflation lambda: %.3f\nMinimum P: %.3g\n", nrow(g), lambda, min(g$P)))
cat(if (lambda < 1) "lambda < 1: possible over-correction\n" else if (lambda <= 1.05) "lambda 1.00-1.05: excellent\n" else if (lambda <= 1.10) "lambda 1.05-1.10: acceptable\n" else "lambda > 1.10: residual stratification - add PCs / check batch and relatedness\n")
EOF_R
Rscript plot_gwas.R | tee "$GWAS_DIR/gwas_summary.txt"
echo "Plots: $GWAS_DIR/plots"
