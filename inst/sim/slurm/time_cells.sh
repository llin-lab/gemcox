#!/bin/bash
# Time replicate 1 of every cell of an experiment (serially, on the current
# node or inside an interactive job) and write slurm/timing_<exp>.tsv:
#   exp  cell  reps  secs_per_replicate  cpus  suggested_minutes
# The suggestion is reps x secs / cpus x 3 + 10 minutes. Writes no results.
#   GEMCOX_CPUS=<cpus per task> slurm/time_cells.sh <exp>
set -euo pipefail
EXP="$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OUT="slurm/timing_${EXP}.tsv"; : > "$OUT"
Rscript run_cell.R --exp "$EXP" --list | while read -r _exp _c cell _rest; do
  Rscript run_cell.R --exp "$EXP" --cell "${cell%:}" --time-only | tail -1 >> "$OUT"
done
column -t "$OUT"
