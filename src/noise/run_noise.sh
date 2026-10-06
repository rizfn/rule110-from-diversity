#!/usr/bin/env bash
# Density of domain walls under noise (Fig. 4C, 4D): 2,048 blocks, 64 runs, 200,000 steps, from a perfect start.
# Output: outputs/noise_NB2048_T200000_seed41/noise_eps<eps>.txt, columns: time, domain walls per block (averaged over the 64 runs)
set -e
cd "$(dirname "$0")"
g++ -O3 -march=native -o walls walls.cpp
NB=2048; T=200000; SEED=41
OUT=outputs/noise_NB${NB}_T${T}_seed${SEED}; mkdir -p $OUT
for e in 0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.003 0.01; do
  ./walls noisec $NB 1 $T $SEED $e > $OUT/noise_eps$e.txt &
done
wait
