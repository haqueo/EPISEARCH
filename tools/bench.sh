#!/usr/bin/env bash
# Time a release build over a fixed dataset at several thread counts.
# Record the numbers before you optimise anything - a speedup you cannot
# quantify is a speedup you cannot defend in an interview.
#
#   ./tools/bench.sh build/release/episearch data.txt 400 61
set -euo pipefail

BIN=${1:?usage: bench.sh <binary> <datafile> <nsamples> <nvars>}
DATA=${2:?}
NSAMPLES=${3:?}
NVARS=${4:?}
REPS=${REPS:-3}

for threads in 1 2 4 8; do
  best=""
  for _ in $(seq "$REPS"); do
    start=$(python3 -c 'import time;print(time.perf_counter())')
    "$BIN" -inputfile "$DATA" -outputfile /dev/null \
           -nsamples "$NSAMPLES" -nvars "$NVARS" \
           -printall n -assoclevel 0.05 -numThread "$threads" >/dev/null 2>&1
    end=$(python3 -c 'import time;print(time.perf_counter())')
    elapsed=$(python3 -c "print(f'{$end-$start:.3f}')")
    if [[ -z "$best" ]] || (( $(python3 -c "print(1 if $elapsed < $best else 0)") )); then
      best=$elapsed
    fi
  done
  printf 'threads=%-2s best of %s: %ss\n' "$threads" "$REPS" "$best"
done
