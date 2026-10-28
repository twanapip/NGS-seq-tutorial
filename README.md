# NGS Variant Calling Pipeline: WGS/WES, Somatic, CNV, Annotation, GWAS

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Shell](https://img.shields.io/badge/shell-bash-blue)
![GATK](https://img.shields.io/badge/GATK-4.x-orange)
![Slurm](https://img.shields.io/badge/HPC-Slurm-lightgrey)

End-to-end, Ready to Run and numbered Bash pipeline for human whole-genome (WGS) and whole-exome (WES)
sequencing: FASTQ → alignment (BWA-MEM) → germline variants (GATK HaplotypeCaller) →
somatic mutations (Mutect2) → copy number (CNVkit) → annotation (SnpEff, VEP, Funcotator,
ANNOVAR) → MAF plots (maftools) → GWAS (PLINK). Runs on a single Linux workstation or a
Slurm cluster. Genome builds hg38 and hg19.

Nothing here is tied to a particular sample: you list your own IDs in text files.

> For research use. Check important findings by an independent method before drawing conclusions.

---

## 1. Folder layout

```
genome-analysis-scripts/
├── README.md
├── LICENSE
├── .gitignore
├── workstation/                 # one Linux machine
│   ├── config.sh                # ALL settings live here
│   ├── lib.sh                   # helper functions (do not edit)
│   ├── run_all.sh               # runs a stage in order
│   └── scripts/00 … 23_*.sh     # numbered steps
└── hpc/                         # SLURM cluster
    ├── config.sh                # same settings + partition/time/mail
    ├── lib.sh
    ├── submit_pipeline.sh       # submits a stage as chained jobs
    └── scripts/                 # same numbered steps + slurm_array/pairs/single.sh
```

Use **one** of `workstation/` or `hpc/`. Always run commands from inside it.

## 2. Quick start

```bash
cd workstation                       # or: cd hpc
nano config.sh                       # set PROJ, BUILD, SEQ_TYPE (and TARGETS for exome)
bash scripts/00_conda_envs.sh        # tools (separate conda envs)
bash scripts/01_reference.sh         # genome + resources + BWA index (hours, once)
bash scripts/03_project_dirs.sh      # folders + empty list files
# WES only:   bash scripts/02_targets.sh
# annotators: bash scripts/11_annotation_db.sh
```

Copy FASTQs into `$PROJ/data/raw/` and fill `$PROJ/lists/samples.txt`, then:

| | Workstation | HPC |
|---|---|---|
| germline calls | `./run_all.sh germline` | `./submit_pipeline.sh germline` |
| joint calling | `./run_all.sh cohort` | `./submit_pipeline.sh cohort` |
| somatic | `./run_all.sh somatic` | `./submit_pipeline.sh somatic` |
| copy number | `./run_all.sh cnv` | `./submit_pipeline.sh cnv` |
| case vs control | `./run_all.sh disease` | `./submit_pipeline.sh disease` |
| MAF + plots | `./run_all.sh maf` | `./submit_pipeline.sh maf` |
| GWAS | `./run_all.sh gwas` | `./submit_pipeline.sh gwas` |
| report | `./run_all.sh review` | `./submit_pipeline.sh review` |

Every script also runs alone and accepts IDs: `bash scripts/05_align.sh <ID> <ID>`. No ID means all IDs in `samples.txt`.
On HPC, scripts 00, 01, 02, 11 and 20 need internet: run them on a login or transfer node.

## 3. Prepare your inputs

**FASTQ names** (in `$PROJ/data/raw/`): `<ID>_R1.fastq.gz` and `<ID>_R2.fastq.gz` (paired-end Illumina).

**Lists** (in `$PROJ/lists/`, one entry per line, `#` = comment). Created empty by script 03.

| File | Content | Used by |
|---|---|---|
| `samples.txt` | every sample ID | all germline steps, joint calling, GWAS |
| `pairs.tsv` | `tumor_id<TAB>normal_id` (`-` for tumor-only). Tumor IDs must be unique | somatic, CNV (tumor), MAF |
| `normals.txt` | normal IDs | pooled normal, custom panel of normals, CNV reference |
| `cases.txt` / `controls.txt` | patient / control IDs | disease-specific step |

IDs must not contain spaces. The ID becomes the BAM sample name (read group `SM`), which Mutect2 matches against.

**Other inputs**: WES capture-kit BED (`TARGETS`), GWAS phenotype file (see section 9).

## 4. What each step does

| # | Script | Does | Main tool | Key output (in `$PROJ`) |
|---|---|---|---|---|
| 00 | `conda_envs` | one conda env per tool group | conda | `ga_core, ga_gatk, ga_snpeff, ga_vep, ga_cnv, ga_gwas, ga_maf` |
| 01 | `reference` | genome, GATK bundle, somatic resources, gene models, BWA index | wget, bwa | `$REF_DIR/<build>/` |
| 02 | `targets` | WES: check BED vs genome, make interval list | GATK | `targets/` |
| 03 | `project_dirs` | folders + list templates | | `lists/` |
| 04 | `fastqc_trim` | read QC, adapter and quality trimming | FastQC, Trim Galore | `qc/fastqc_raw`, `data/trimmed` |
| 05 | `align` | map reads, sorted BAM | BWA-MEM, samtools | `bam/<ID>.sorted.bam` |
| 06 | `markdup_bqsr` | mark duplicates, recalibrate base qualities | GATK | `bam/<ID>.recal.bam` |
| 07 | `alignment_qc` | mapping, insert size, coverage (WGS metrics or WES hybrid-selection) | GATK, samtools | `qc/metrics/` |
| 08 | `call_germline` | per-sample GVCF, then VCF | HaplotypeCaller | `gvcf/`, `vcf/<ID>.raw.vcf.gz` |
| 09 | `filter_variants` | hard filters (or VQSR) | GATK | `vcf/<ID>.pass.vcf.gz` |
| 10 | `joint_calling` | all samples together, then filter | GenomicsDB | `cohort/cohort.pass.vcf.gz` |
| 11 | `annotation_db` | download databases for chosen annotators | | |
| 12 | `annotate` | add genes, effects, frequencies, clinical info | SnpEff, Funcotator, VEP, ANNOVAR | `annotation/<vcf>/` |
| 13 | `variant_stats_qc` | stats, one-page QC summary, table, BED, coverage track | bcftools | `qc/<ID>.qc_summary.txt` |
| 14 | `somatic_prepare` | pooled normal BAM or custom panel of normals | samtools, GATK | `somatic/<mode>/` |
| 15 | `somatic_mutect2` | call somatic mutations | Mutect2 | `somatic/<mode>/<tumor>.raw.vcf.gz` |
| 16 | `somatic_filter` | contamination, orientation bias, filters, table | GATK, bcftools | `*.high_confidence.vcf.gz`, `*.mutations.tsv` |
| 17 | `cnv` | copy number calls, filters, exports, plots | CNVkit | `cnv/<mode>/` |
| 18 | `disease_variants` | patient-specific and case-enriched variants | bcftools | `cohort/disease/` |
| 19 | `maf_visualize` | VCF to MAF, oncoplot, burden, signatures, drivers | vcf2maf, maftools | `maf/plots/` |
| 20 | `gwas_prepare` | VCF to PLINK files (or 1000 Genomes practice data) | bcftools, PLINK | `gwas/data/` |
| 21 | `gwas_qc_pca` | sample/variant QC, ancestry PCA | PLINK, R | `gwas/results/clean.*`, `pca.*` |
| 22 | `gwas_assoc_plot` | association test, Manhattan, QQ, lambda | PLINK, R | `gwas/plots/` |
| 23 | `review_outputs` | MultiQC and file counts | MultiQC | `qc/multiqc/` |

## 5. Concepts in one minute

- **FASTQ** = raw reads. **BAM** = reads placed on the genome. **VCF** = list of differences from the reference. **GVCF** = VCF that also records "no variant here", so samples can be combined later.
- **WGS vs WES**: WGS reads the whole genome (~30x, finds all variant types, good for CNVs). WES reads only exons via a capture kit (~100x, cheaper, coding variants only, weaker CNVs). Always use the *exact* BED of your kit.
- **Trimming** removes adapters and poor ends so reads align cleanly.
- **Duplicates** are copies made by PCR; marking them stops one molecule being counted many times.
- **BQSR** corrects systematic over/under-confident base qualities using known variant sites.
- **Hard filter vs VQSR**: hard filters apply fixed cutoffs and work for any sample count. VQSR learns from truth sets and needs many samples (roughly 30+), so it is only offered for the joint cohort.
- **Germline vs somatic**: germline variants are inherited and in every cell (about 50% or 100% of reads). Somatic variants arise later (tumor), can be as low as 5-30% of reads, and need a normal for comparison.
- **Somatic strategies** (set `SOMATIC_MODE`):

| Mode | When | Specificity | Notes |
|---|---|---|---|
| `matched` | tumor + normal from the same person | highest | gold standard |
| `pooled` | several unrelated normals, no matched one | high | normals merged into one reference |
| `custom_pon` | 40+ normals from your own kit/platform | good | removes your own technical artifacts |
| `tumor_only` | no normal at all | lowest | germline leaks in, strict filters, validate everything |

- **Panel of Normals (PoN)** = sites that look "mutated" in many normals, so they are artifacts. **gnomAD AF** = population frequency, used to spot germline variants. **Contamination** = DNA from another person (above 5% hurts sensitivity). **Orientation bias** = chemical artifacts on one strand (common in FFPE).
- **CNV** = a segment of DNA gained or lost. Measured as log2 ratio of read depth vs a reference (0 = 2 copies, -1 = 1 copy, +0.58 = 3 copies). Needs a good normal reference to cancel technical bias.
- **Annotation** turns "chr17:43,044,295 G>A" into gene, protein change, consequence, frequency and clinical significance. Impact classes: HIGH (stop, frameshift, splice), MODERATE (missense), LOW (synonymous), MODIFIER (non-coding).
- **GWAS** tests every common variant for association with a trait. **Population stratification** (ancestry differences) causes false hits, so we add principal components (PCs) as covariates. **Lambda** (genomic inflation) near 1.0 means well controlled. Genome-wide significance is p < 5e-8.

## 6. Settings that matter (`config.sh`)

| Setting | Values | Meaning |
|---|---|---|
| `BUILD` | `hg38`, `hg19` | genome build (hg19 = Broad b37, contigs `1..22,X,Y`) |
| `SEQ_TYPE` | `WGS`, `WES` | WES also needs `TARGETS` (kit BED) and adds depth filters |
| `FILTER_MODE` | `hard`, `vqsr` | hard = default; vqsr only for large cohorts |
| `ANNOTATORS` | any of `snpeff funcotator vep annovar` | tools used by script 12 (default only snpeff) |
| `SOMATIC_MODE` | `matched pooled custom_pon tumor_only` | see table above |
| `CNV_MODE` | `tumor`, `germline` | reference and thresholds |
| `TRAIT_TYPE` | `quant`, `binary` | linear vs logistic regression |
| `N_PCS` | 5-10 typical | PCs used as GWAS covariates (use the scree plot elbow) |
| `CLEAN_INTERMEDIATE` | 0/1 | delete sorted/dedup BAMs to save disk |
| `CORES`, `MEM_GB` | numbers | threads and Java RAM |

Annotation tool choice: **SnpEff** is fast and simple; **VEP** is the most detailed consequence predictor and can add SIFT/PolyPhen/CADD via `DBNSFP_FILE`; **Funcotator** keeps everything inside GATK and has germline and somatic data sources; **ANNOVAR** is strong for population frequencies (free registration needed, set `ANNOVAR_DIR`).

## 7. Thresholds used (all editable in the scripts or config)

**Germline hard filters** (SNPs and indels filtered separately)

| | QD | FS | MQ | SOR | MQRankSum | ReadPosRankSum |
|---|---|---|---|---|---|---|
| SNP | < 2 | > 60 | < 40 | > 3 | < -12.5 | < -8 |
| Indel | < 2 | > 200 | n/a | > 10 | n/a | < -20 |

WES adds `DP < WES_MIN_DP` and `DP > WES_MAX_DP`. VQSR sensitivity: SNP 99.0, indel 95.0.

**Somatic high-confidence filter** (after `FilterMutectCalls` PASS)

| Mode | Tumor AF | Tumor DP | TLOD | Normal |
|---|---|---|---|---|
| matched | >= 0.05 | >= 10 | >= 6.3 | DP >= 10 and AF <= 0.03 |
| pooled | >= 0.06 | >= 12 | >= 6.3 | DP >= 10 and AF <= 0.03 |
| custom_pon | >= 0.05 | >= 10 | >= 6.3 | none |
| tumor_only | >= 0.10 | >= 20 | >= 10 | population AF <= 0.001 (`POPAF >= 3`) |

**CNV** call thresholds on log2: germline `-1.1,-0.4,0.3,0.7`; tumor `-1.0,-0.5,0.5,1.0`. Kept if copy number is not 2, size >= 10 kb and >= 3 bins. Germline calls overlapping common CNVs (DGV, 50%) are removed to `*.rare.bed`.

**Disease-specific**: control AF < 0.05, gnomAD AF < 0.01, HIGH or MODERATE impact; case-enriched also needs case AF > control AF and seen in >= 2 cases.

**GWAS QC**: sample missingness < 5%, variant missingness < 5%, MAF >= 0.01, HWE p > 1e-6. Optional: heterozygosity within 3 SD, remove relatives.

## 8. QC benchmarks (what good looks like)

| Metric | Good |
|---|---|
| Reads aligned | > 95% |
| Duplicate rate | < 30% (WGS; lower is better) |
| Mean coverage | WGS >= 30x; WES >= 100x clinical, >= 50x research |
| WES reads on target | > 80%; targets at >= 20x > 95% |
| Base quality / adapters | Q > 20 across reads; adapters < 5% |
| SNPs per genome (WGS) | 4-5 million (exome 20-25 thousand) |
| Indels (WGS) | 0.5-0.8 million |
| Ti/Tv | ~2.0-2.1 (WGS), ~2.8-3.3 (exome); somatic 2-3 |
| Het/Hom | 1.5-2.0 |
| In dbSNP | > 95% |
| Contamination (somatic) | < 2% excellent, 2-5% fine, > 5% poor |
| Somatic counts | matched 500-5,000 (WGS); < 100 or > 10,000 is a red flag |
| GWAS lambda | 1.00-1.05 excellent, up to 1.10 acceptable |
| CNV | smooth baseline near 0; sex chromosomes match expected sex |

Typical hardware: 64 GB RAM minimum for human WGS (128 GB better), ~100 GB FASTQ and ~100 GB BAM per 30x genome, SSD helps, multi-core. Windows or laptops are not suitable for WGS/WES; use Linux or a cluster.

## 9. GWAS notes

- Real GWAS needs hundreds to thousands of people. Small cohorts only teach the workflow.
- **Own data**: set `GWAS_VCF` to a joint-called cohort VCF. **Practice**: leave it empty with `GWAS_DATA="1000g"` (downloads 2,504 genomes, GRCh37). `SIMULATE_PHENO=1` makes a fake trait (ancestry + sex + noise, no true genetics), so any hit is chance.
- `PHENO_FILE`: no header, three columns `FID IID value` (binary: 1 control, 2 case, -9 missing). `COVAR_EXTRA`: file *with* header (`FID IID age batch ...`); sex and PCs are added automatically.
- Replication in an independent cohort and meta-analysis (e.g. METAL: combine effect size, SE and sample size per variant, fixed effects if I² is low) are needed before believing a hit. For a single gene region, subset with PLINK `--chr --from-bp --to-bp`. PRS methods use many sub-threshold variants, not only significant ones.

## 10. Aligner chooser (this repo uses BWA-MEM)

| Data | Use |
|---|---|
| Short-read DNA, variant calling | **BWA-MEM** |
| ChIP-seq / ATAC-seq | Bowtie2 |
| Bulk RNA-seq, need coordinates/junctions | STAR (30+ GB RAM) or HISAT2 (4-8 GB) |
| RNA-seq gene counts only | Salmon or Kallisto |
| 10x single-cell | STARsolo, Cell Ranger, Alevin-fry, kb-python |
| Long reads (PacBio/Nanopore) | Minimap2 (`map-ont`, `map-pb`, or `splice` for RNA); NGMLR for structural variants |
| Reference is a transcriptome FASTA | a DNA aligner is fine (no introns) |

## 11. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `conda not found` | conda not on PATH in scripts | set `CONDA_BASE` (or `HPC_INIT` on HPC) |
| GATK: missing read group | aligned without `-R` | re-run script 05 |
| Mutect2: normal sample not found | ID differs from BAM sample name | IDs in lists must equal the IDs used for alignment |
| Low mapping rate | wrong build/species, adapters | check `BUILD`, FastQC, `samtools flagstat` |
| WES: contig error in step 02 | `chr1` vs `1` | use a BED matching your build |
| Java version clash | tools sharing an env | scripts already use separate envs; do not merge them |
| VQSR fails | too few variants/samples | `FILTER_MODE="hard"` |
| Out of memory | heap too large or too small | adjust `MEM_GB`; keep below machine RAM |
| Many somatic calls, few PASS | normal and artifacts | expected: raw calls are ~1-2% of final |
| Too few GWAS hits / inflated lambda | power / stratification | more samples; more PCs; relatedness and batch checks |
| Download failed (hg19 file) | bucket file name changed | check the bucket listing and edit `lib.sh` |
| Disk full | BAMs are big | `CLEAN_INTERMEDIATE=1`, delete `tmp/` |

Speed tips: WGS HaplotypeCaller can be run per chromosome with `-L` and merged with `MergeVcfs`; joint calling and PoN building already work per chromosome. Scatter Mutect2 the same way for exomes (merge VCFs with `MergeVcfs`, stats with `MergeMutectStats`, filter once).

## 12. Design Notes (bugs avoided)

- Normal-sample columns are found **by name**, never assumed (`AF[1:0]` = normal, not `AF[0:1]`), and normal depth is required.
- Tumor-only germline filter uses `POPAF >= 3` (POPAF is -log10 of frequency), not `< 0.001`.
- A pooled normal BAM gets one read-group sample name, otherwise Mutect2 cannot find it.
- Separate SNP and indel hard filters; VQSR is not used on single samples or exomes.
- GATK4 argument styles (`--BAIT_INTERVALS`), correct Funcotator downloader flags, no comments after line-continuation backslashes, `$HOME` instead of `~` inside quotes.
- WES: BQSR learns on targets, but the final BAM keeps off-target reads for CNV.
- CNV filters pick columns by header (`cn`, `probes`), not by position.
- MAF signature step uses the BSgenome matching your build.
- Frequency tags are added from lookup tables (not tag renaming), one ALT per line.
- Each tool group has its own conda env (GATK and SnpEff need different Java versions).

## 13. Credits and license

Built on GATK, BWA, samtools/bcftools, FastQC, Trim Galore, SnpEff, Ensembl VEP, ANNOVAR, CNVkit, vcf2maf, maftools, PLINK, and the Broad/1000 Genomes resources. Cite the tools you use. MIT license, see `LICENSE`.
