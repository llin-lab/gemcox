#!/bin/bash
# Submit every MISSING cell of an experiment as its own one-task array job
# (array index = cell id), each with its own time limit. Completed cells are
# skipped, so rerunning this script is safe.
#
#   slurm/submit.sh <exp> [cpus] [partition]
#
# Time limits come from slurm/timing_<exp>.tsv, made by slurm/time_cells.sh
# (one replicate per cell). Without it the default is 24:00:00 per cell.
# Run from gemcox/inst/sim. Logs: slurm/logs/<exp>_cell<k>_<jobid>.log
set -euo pipefail
EXP="$1"; CPUS="${2:-16}"; PART="${3:-}"
SIM_DIR="$(pwd)"
[ -f run_cell.R ] || { echo "run from gemcox/inst/sim"; exit 1; }
mkdir -p slurm/logs
TIMING="slurm/timing_${EXP}.tsv"
Rscript run_cell.R --exp "$EXP" --list | while read -r _exp _c cell _reps _r status; do
  cell="${cell%:}"
  if [ "$status" = "complete" ]; then echo "cell $cell complete: skipped"; continue; fi
  mins=1440
  if [ -f "$TIMING" ]; then
    m=$(awk -v c="$cell" '$2 == c {print $6}' "$TIMING")
    [ -n "$m" ] && mins="$m"
  fi
  printf -v tlim "%02d:%02d:00" $((mins / 60)) $((mins % 60))
  extra=(); [ -n "$PART" ] && extra+=(--partition="$PART")
  jid=$(sbatch --parsable --job-name="gemcox_${EXP}" --array="$cell" --cpus-per-task="$CPUS" \
        --time="$tlim" --mem-per-cpu=2G --output="slurm/logs/${EXP}_cell%a_%A.log" \
        --export=ALL,GEMCOX_SIM_DIR="$SIM_DIR" "${extra[@]}" slurm/run_task.sh "$EXP")
  echo "cell $cell submitted: job $jid, time $tlim, $CPUS cpus"
done
