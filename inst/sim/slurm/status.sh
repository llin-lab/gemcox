#!/bin/bash
# Status of every cell of an experiment: complete (cell file present),
# running or pending (in squeue), failed (a log ends in an error without a
# cell file), or missing.   slurm/status.sh <exp>     (from gemcox/inst/sim)
set -uo pipefail
EXP="$1"
Q=""; command -v squeue >/dev/null && Q=$(squeue -h -u "$USER" -n "gemcox_${EXP}" -r -o "%K %T" 2>/dev/null || true)
Rscript run_cell.R --exp "$EXP" --list | while read -r _exp _c cell _reps _r status; do
  cell="${cell%:}"
  if [ "$status" = "complete" ]; then echo "cell $cell: complete"; continue; fi
  st=$(echo "$Q" | awk -v c="$cell" '$1 == c {print tolower($2)}' | head -1)
  if [ -n "$st" ]; then echo "cell $cell: $st"; continue; fi
  log=$(ls -t slurm/logs/${EXP}_cell${cell}_*.log 2>/dev/null | head -1 || true)
  if [ -n "$log" ] && grep -qE "Error|Execution halted|CANCELLED|TIME LIMIT|oom" "$log"; then
    echo "cell $cell: FAILED (see $log)"
  else
    echo "cell $cell: missing"
  fi
done
