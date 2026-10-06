# Code for "A repeating sequence of simple rules creates a universal computer"

Everything needed to regenerate the figures and the numbers in the paper.

## Requirements
- Python 3 with NumPy and Matplotlib (`uv sync` from the repository root)
- gcc and g++ (only to regenerate the simulation data)
- Optional: the Liberation Sans font (the figures fall back to Arial or the default sans otherwise)

## Layout

Each folder holds the code for one part of the paper. Simulation data goes to `<folder>/outputs/` (not tracked by git),
figures to `<folder>/plots/`. C/C++ codes are compiled next to their source. Every script can be run from anywhere.

```
common/                     shared Python modules
  mixture.py                the 36|37|5 mixture, the block code, elementary rules
  beats.py                  beat (phase) of each block for the colour maps; hand-made starts for Fig. 3B
  style.py                  figure fonts and colours
construction/               Sec. 2, Figs. 1-2
  verify_construction.py    checks that 3 mixture steps = 1 rule-110 step for every configuration
  fig1.py, fig2.py          -> plots/
random_starts/              Sec. 3, Fig. 3
  lockcharge.c              lock times from random starts, with the walls' net charge
  run_lockcharge.sh         -> outputs/lock_charge/charge_N<N>_tmax<TMAX>_s<seed>.txt  ("charge lock_time" per run)
  lock_stats.py             lock fractions, median lock times and scaling exponents
  fig3.py                   -> plots/
  outputs/lock_times/       lock times, N = 32-4096 (one value per run; -1 = not locked), from an earlier code
noise/                      Sec. 4, Fig. 4
  walls.cpp                 walls created by single flips and by noise
  run_noise.sh              -> outputs/noise_NB2048_T200000_seed41/noise_eps<eps>.txt  (time, walls per block)
  measure_parameters.py     escape probability p and wall speeds v_R, v_L; results are stored in params.py
  params.py                 the measured p, v_R, v_L used in Fig. 4
  fig4.py                   -> plots/
info_cost/                  Discussion
  info_cost.py              description lengths of the rules and the mixture
```

## Reproducing the paper

From existing data (a few minutes):
```
python src/construction/verify_construction.py
python src/construction/fig1.py; python src/construction/fig2.py
python src/random_starts/fig3.py; python src/random_starts/lock_stats.py
python src/noise/fig4.py
python src/info_cost/info_cost.py
```

Regenerating the data:
```
src/noise/run_noise.sh                      # minutes
src/random_starts/run_lockcharge.sh         # hours on many cores; the 16,384-block runs take several hours per job
python src/noise/measure_parameters.py      # about a minute
```

## Notes
- The C codes simulate 64 independent lattices at once, one per bit of a 64-bit word.
- `random_starts/outputs/lock_times` (Fig. 3C1, N <= 4096) was produced with an earlier bit-parallel code with the same
  dynamics and cannot be regenerated exactly. `run_lockcharge.sh` regenerates statistically equivalent data (its second
  column is the lock time) with different random seeds. The 8,192- and 16,384-block points in Fig. 3C1 come from
  `outputs/lock_charge`.
- The figures in `plots/` are identical to those in the paper.
