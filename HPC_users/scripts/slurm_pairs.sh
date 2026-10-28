#!/bin/bash
# Generic array job: task N runs STEP for the Nth tumor in pairs.tsv.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
[ -z "$HPC_INIT" ] || eval "$HPC_INIT"
T="$(tumors | sed -n "${SLURM_ARRAY_TASK_ID}p")"
[ -n "$T" ] || { echo "No tumor for array task $SLURM_ARRAY_TASK_ID"; exit 1; }
bash "$SLURM_SUBMIT_DIR/$STEP" "$T" "$@"
