#!/usr/bin/env bash
# High-confidence somatic VCFs -> MAF (vcf2maf + VEP) -> merged MAF -> maftools plots (oncoplot, lollipop, burden, Ti/Tv,
# signatures, drivers, domains, interactions).  Output: maf/ and maf/plots/
# Usage: bash scripts/19_maf_visualize.sh [tumor_id ...]      (uses SOMATIC_MODE; VEP cache from script 11 required)
set -euo pipefail
source "$(dirname "$0")/../config.sh"
OUT="$PROJ/maf"; mkdir -p "$OUT/per_sample" "$OUT/plots"
V2M="$REF_DIR/tools/vcf2maf"
if [ ! -d "$V2M" ]; then mkdir -p "$REF_DIR/tools"; git clone https://github.com/mskcc/vcf2maf.git "$V2M"; fi
use_env vep
VEP_PATH="$(dirname "$(command -v vep)")"
for T in $(tumors_from_args "$@"); do
  N="$(normal_of "$T")"; NIDS=()
  case "$SOMATIC_MODE" in
    matched) NIDS=(--normal-id "$N" --vcf-normal-id "$N") ;;
    pooled)  NIDS=(--normal-id pooled_normal --vcf-normal-id pooled_normal) ;;
  esac
  zcat "$PROJ/somatic/$SOMATIC_MODE/$T.high_confidence.vcf.gz" > "$OUT/per_sample/$T.vcf"
  perl "$V2M/vcf2maf.pl" --input-vcf "$OUT/per_sample/$T.vcf" --output-maf "$OUT/per_sample/$T.maf" \
    --tumor-id "$T" --vcf-tumor-id "$T" "${NIDS[@]}" --ncbi-build "$VEP_ASSEMBLY" --species homo_sapiens \
    --ref-fasta "$REF_FASTA" --vep-path "$VEP_PATH" --vep-data "$VEP_CACHE" --cache-version "$VEP_CACHE_VERSION"
  rm -f "$OUT/per_sample/$T.vcf"
done

if [ "$SEQ_TYPE" = "WES" ]; then CALLABLE_MB=$(awk '{s+=$3-$2} END{printf "%.1f", s/1e6}' "$TARGETS"); else CALLABLE_MB=3088; fi
export MAF_DIR="$OUT" CALLABLE_MB BSGENOME
use_env maf
cat > "$OUT/plot_maf.R" <<'EOF_R'
suppressPackageStartupMessages({ library(maftools); library(ggplot2); library(dplyr) })
out <- Sys.getenv("MAF_DIR"); mb <- as.numeric(Sys.getenv("CALLABLE_MB")); genome <- Sys.getenv("BSGENOME")
setwd(out)
files <- Sys.glob(file.path(out, "per_sample", "*.maf"))
if (length(files) == 0) stop("No per-sample MAF files found")
maf <- if (length(files) > 1) merge_mafs(files, verbose = TRUE) else read.maf(files[1])
write.table(maf@data, "merged_mutations.maf", sep = "\t", quote = FALSE, row.names = FALSE)

safe <- function(file, w, h, expr) {
  tryCatch({ pdf(file.path("plots", file), width = w, height = h); force(expr); dev.off() },
           error = function(e) { while (dev.cur() > 1) dev.off(); message("SKIPPED ", file, ": ", conditionMessage(e)) })
}
safe("summary.pdf", 12, 8, plotmafSummary(maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE))
safe("oncoplot_top20.pdf", 12, 8, oncoplot(maf, top = 20, removeNonMutated = FALSE))
genes <- head(getGeneSummary(maf)$Hugo_Symbol, 5)
safe("lollipop_top5.pdf", 12, 8, for (g in genes) lollipopPlot(maf, gene = g, showMutationRate = TRUE))

ss <- getSampleSummary(maf) %>% mutate(mut_per_mb = total / mb)
write.table(ss, "sample_mutation_burden.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
safe("mutation_burden.pdf", 10, 6, print(
  ggplot(ss, aes(reorder(Tumor_Sample_Barcode, -mut_per_mb), mut_per_mb)) + geom_col(fill = "#1F78B4") +
    geom_hline(yintercept = mean(ss$mut_per_mb), linetype = "dashed", colour = "red") +
    labs(x = "Sample", y = "Mutations per Mb", title = "Mutation burden") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))))
gs <- head(getGeneSummary(maf), 20)
safe("top_mutated_genes.pdf", 8, 8, print(
  ggplot(gs, aes(reorder(Hugo_Symbol, MutatedSamples), MutatedSamples)) + geom_col(fill = "#33A02C") + coord_flip() +
    labs(x = "Gene", y = "Mutated samples", title = "Most frequently mutated genes") + theme_minimal()))

safe("titv.pdf", 10, 6, plotTiTv(res = titv(maf, plot = FALSE, useSyn = TRUE)))
safe("apobec_signature.pdf", 12, 8, {
  has_chr <- any(grepl("^chr", maf@data$Chromosome))
  tnm <- trinucleotideMatrix(maf, prefix = if (has_chr) NULL else "chr", add = TRUE, ref_genome = genome)
  plotApobecDiff(tnm = tnm, maf = maf, pVal = 0.2)
})
safe("driver_genes.pdf", 10, 8, { od <- oncodrive(maf, AACol = "HGVSp_Short", minMut = 5); plotOncodrive(od, fdrCutOff = 0.1, useFraction = TRUE) })
safe("protein_domains.pdf", 12, 8, pfamDomains(maf, top = 10))
safe("somatic_interactions.pdf", 12, 10, somaticInteractions(maf, top = 25, pvalue = c(0.05, 0.1)))
message("Done: ", file.path(out, "plots"))
EOF_R
Rscript "$OUT/plot_maf.R"
