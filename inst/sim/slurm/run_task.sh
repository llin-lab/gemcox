#!/bin/bash
# One SLURM array task = one cell. Submitted by slurm/submit.sh; the array
# index IS the cell id. Usage (via sbatch): run_task.sh <exp>
set -euo pipefail
EXP="$1"
CELL="${SLURM_ARRAY_TASK_ID:?run this through sbatch --array}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export GEMCOX_WORKERS="${SLURM_CPUS_PER_TASK:-1}"
cd "${GEMCOX_SIM_DIR:?set GEMCOX_SIM_DIR to gemcox/inst/sim (submit.sh does)}"
echo "[$(date)] $EXP cell $CELL on $(hostname), $GEMCOX_WORKERS workers, commit $(git rev-parse --short HEAD)"
Rscript run_cell.R --exp "$EXP" --cell "$CELL"
echo "[$(date)] $EXP cell $CELL finished"
