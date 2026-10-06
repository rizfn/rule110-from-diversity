#!/usr/bin/env bash
# Lock times from random starts (Fig. 3C1, 3C2 and the lock fractions in Sec. 3).
# Each job = 64 runs; 8 seeds per size = 512 runs. Jobs run in parallel on all cores (CORES=8 ./run_lockcharge.sh to limit).
# Output: outputs/lock_charge/charge_N<N>_tmax<TMAX>_s<seed>.txt, one line per run: "charge lock_time"
#   charge 0 = balanced, < 0 = excess left-moving walls, > 0 = excess right-moving walls; lock_time -1 = not locked
# Note: the paper's C1 points for N <= 4096 (outputs/lock_times) were produced by an earlier bit-parallel code with the same
# dynamics; this script regenerates statistically equivalent data. The largest sizes take hours per job.
set -e
cd "$(dirname "$0")"
gcc -O3 -march=native -o lockcharge lockcharge.c
OUT=outputs/lock_charge; mkdir -p $OUT
CORES=${CORES:-$(nproc)}
jobs=()
#        N      TMAX (steps)
for spec in "128 2000000" "256 2000000" "512 4000000" "1024 8000000" "2048 25000000" \
            "4096 60000000" "8192 150000000" "16384 400000000"; do
  set -- $spec
  for s in 1 2 3 4 5 6 7 8; do jobs+=("$1 $2 $s"); done
done
printf '%s\n' "${jobs[@]}" | xargs -P "$CORES" -L 1 sh -c "./lockcharge \$0 \$1 \$2 > $OUT/charge_N\$0_tmax\$1_s\$2.txt"
echo "done: $(ls $OUT/charge_N*_s*.txt | wc -l) files"
