#!/bin/bash
# Generic single job: runs STEP once (cohort, CNV, GWAS, ...).
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
[ -z "$HPC_INIT" ] || eval "$HPC_INIT"
bash "$SLURM_SUBMIT_DIR/$STEP" "$@"
