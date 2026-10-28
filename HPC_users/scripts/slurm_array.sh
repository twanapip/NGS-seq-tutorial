#!/bin/bash
# Generic array job: task N runs STEP for the Nth ID in samples.txt.  Submitted by submit_pipeline.sh.
set -euo pipefail
source "$SLURM_SUBMIT_DIR/config.sh"
[ -z "$HPC_INIT" ] || eval "$HPC_INIT"
S="$(read_list "$SAMPLE_LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")"
[ -n "$S" ] || { echo "No sample for array task $SLURM_ARRAY_TASK_ID"; exit 1; }
bash "$SLURM_SUBMIT_DIR/$STEP" "$S" "$@"
