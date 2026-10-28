# HPC (SLURM) config. Edit, then run from the hpc/ directory: setup scripts on a login node, the rest via ./submit_pipeline.sh
# Sample lists live in $PROJ/lists (created by scripts/03_project_dirs.sh) - see README "Prepare your inputs".

# ---------------- paths ----------------
export PROJ="$HOME/genome_analysis"      # all results go here (needs lots of disk space)
export REF_DIR="$HOME/references"        # genome + resource files (can be shared by projects)
export ENV_PREFIX="ga"                    # conda envs are named ga_core, ga_gatk, ...
export CONDA_BASE=""                      # leave empty to auto-detect conda

# ---------------- data type ----------------
export BUILD="hg38"                       # hg38, or hg19 (= Broad "b37": contigs 1..22,X,Y without "chr")
export SEQ_TYPE="WGS"                     # WGS or WES
export TARGETS=""                         # WES only: capture-kit BED from your kit vendor. WGS: leave empty
export PAD=100                            # WES only: bases added around each target when calling

# ---------------- compute ----------------
export CORES="${SLURM_CPUS_PER_TASK:-16}"      # per-job CPUs (set by sbatch --cpus-per-task)
export MEM_GB=64                          # Java RAM per job; each job requests MEM_GB+8 GB
export CLEAN_INTERMEDIATE=0               # 1 = delete sorted/dedup BAMs after the final recalibrated BAM exists

# ---------------- input lists (one ID per line, # = comment) ----------------
export SAMPLE_LIST="$PROJ/lists/samples.txt"      # every sample ID (FASTQs: data/raw/<ID>_R1.fastq.gz, <ID>_R2.fastq.gz)
export PAIRS_FILE="$PROJ/lists/pairs.tsv"         # somatic: tumor_id <TAB> normal_id   (normal "-" = tumor-only)
export NORMALS_FILE="$PROJ/lists/normals.txt"     # normal IDs for pooled normal / custom PoN / CNV reference
export CASES_FILE="$PROJ/lists/cases.txt"         # disease-specific step: patient IDs
export CONTROLS_FILE="$PROJ/lists/controls.txt"   # disease-specific step: control IDs
export USE_TRIMMED=1                              # 1 = align trimmed reads, 0 = align raw reads

# ---------------- germline calling / filtering ----------------
export FILTER_MODE="hard"                 # hard = safest for single samples, small cohorts and exomes; vqsr = large cohorts
export VQSR_SNP_SENS=99.0                 # VQSR truth-sensitivity levels (only used when FILTER_MODE=vqsr)
export VQSR_INDEL_SENS=95.0
export WES_MIN_DP=20                      # WES only: extra depth filters
export WES_MAX_DP=500

# ---------------- annotation ----------------
ANNOTATORS=(snpeff)                       # any of: snpeff funcotator vep annovar
export SNPEFF_DB=""                       # empty = GRCh38.105 (hg38) or GRCh37.75 (hg19)
export SNPEFF_MEM=16
export VEP_CACHE_VERSION=115
export VEP_CACHE="$REF_DIR/vep_cache"
export FUNCO_DIR="$REF_DIR/funcotator"
export ANNOVAR_DIR=""                     # folder with annotate_variation.pl (needs free registration, see README)
export ANNOVAR_PROTOCOLS="refGene,gnomad312_genome,1000g2015aug_all,clinvar_20220320,dbnsfp42c"
export ANNOVAR_OPERATIONS="g,f,f,f,f"
export DBNSFP_FILE=""                     # optional: dbNSFP .gz file to add SIFT/PolyPhen/CADD through VEP

# ---------------- somatic ----------------
export SOMATIC_MODE="matched"             # matched | pooled | custom_pon | tumor_only
export MAKE_BEDGRAPH=1                    # 1 = write coverage bedGraph in the QC step (large file for WGS)

# ---------------- copy number (CNVkit) ----------------
export CNV_MODE="tumor"                   # tumor (vs normals.txt reference) or germline
export CNV_MIN_SIZE=10000                 # keep calls >= this many bp
export CNV_MIN_PROBES=3                   # keep calls supported by >= this many bins
export CNV_THRESH_GERMLINE="-1.1,-0.4,0.3,0.7"
export CNV_THRESH_TUMOR="-1.0,-0.5,0.5,1.0"
export PURITY=""                          # tumor purity 0-1 if known (empty = not used)
export DGV_URL="http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt"   # common germline CNVs (hg38)
export CNV_CLINVAR=0                      # 1 = overlap tumor CNVs with ClinVar CNV records

# ---------------- disease / patient-specific variants ----------------
export DISEASE_CONTROL_AF_MAX=0.05
export DISEASE_GNOMAD_AF_MAX=0.01
export DISEASE_MIN_CASE_AC=2

# ---------------- GWAS (PLINK) ----------------
export GWAS_VCF=""                        # your joint-called cohort VCF (.vcf.gz). Empty + GWAS_DATA=1000g = practice data
export GWAS_DATA="1000g"                  # 1000g (downloads 2,504 genomes, GRCh37) or none
export PHENO_FILE=""                      # text file, NO header: FID IID phenotype (binary: 1=control 2=case, -9 missing)
export TRAIT_TYPE="quant"                 # quant (linear) or binary (logistic)
export COVAR_EXTRA=""                     # optional file WITH header: FID IID age batch ... (sex + PCs are added automatically)
export N_PCS=10
export SIMULATE_PHENO=0                   # 1 = make a practice trait (no real genetics) when you have no phenotype
export GWAS_MIND=0.05                     # drop samples with > 5% missing genotypes
export GWAS_GENO=0.05                     # drop variants with > 5% missing
export GWAS_MAF=0.01
export GWAS_HWE=1e-6
export HET_FILTER=0                       # 1 = drop samples with heterozygosity > 3 SD from the mean
export RELATEDNESS_FILTER=0               # 1 = drop one of each related pair (recommended for real cohorts)
export POP_PANEL=""                       # optional: file with columns sample, super_pop to colour the PCA plot

# ---------------- SLURM ----------------
export PARTITION="normal"                 # your cluster's partition
export TIME="2-00:00:00"                  # wall time per job
export MAIL="you@example.com"             # failure notices
export HPC_INIT=""                        # e.g. "module load miniconda3" (run before conda is used)
export CPUS_LIGHT=4                       # CPUs for light steps (QC summary, filters, annotation)

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
