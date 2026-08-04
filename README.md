# GENESIS GPU Benchmark

This repository contains input templates and Python scripts for benchmarking
the GENESIS `spdyn` GPU molecular-dynamics implementation.

The benchmark runs a matrix of systems, ensembles, and time steps. For each
selected cell it can:

1. run a GENESIS autotuning pass,
2. write a pinned production input,
3. run warmups,
4. run repeated production measurements,
5. save raw logs and a detailed CSV file.

## Requirements

- Python 3.10 or newer
- `mpirun`
- a GENESIS `spdyn` executable
- the compressed input archives in `data/*.tgz`

The command examples below assume you run them from this repository root.

## Quick Start

Run one small benchmark:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --ensembles nve \
  --dt 2 \
  --warmup 1 \
  --measure 3
```

Run the default benchmark matrix:

```bash
python run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn
```

The default run uses:

- all systems: `dhfr,apoa1,uun,factorix,bpti,dppc,ake,stmv,cellulose`
- all ensembles: `nve,nvt,npt`
- both time steps: `2fs,4fs` for every system
- kernel autotuning
- `100000` production steps
- `50000` autotuning steps
- `eneout_period = 1000`
- `1` MPI process
- `1` OpenMP thread

## Common Examples

Run the 23,558-atom AMBER JAC/DHFR NVE benchmark at both time steps:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --ensembles nve \
  --dt 2,4
```

Run DHFR and cellulose:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr,cellulose \
  --ensembles nve \
  --dt 2
```

Use shorter runs for a smoke test:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --ensembles nve \
  --dt 2 \
  --warmup 0 \
  --measure 1 \
  --nsteps 1000 \
  --tune-nsteps 1000
```

Set MPI and OpenMP parallelism:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --mpi-procs 1 \
  --omp-threads 1
```

Override the energy-output period:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --nsteps 20000 \
  --tune-nsteps 10000 \
  --eneout-period 1000
```

Both `--nsteps` and `--tune-nsteps` must be exact multiples of
`--eneout-period`.

## Output

Each benchmark creates:

```text
results/<timestamp>.csv
results/<timestamp>/
  benchmark.log
  summary.log
  inputs/
  autotune/
  production/
```

`benchmark.log` contains the full progress log.

`summary.log` contains only the final table:

```text
=== ns/day (mean/median +- std, cv%) : tuners=kernel, measure=10, nsteps=100000 ===
system     ens   dt           mean       median      +-std     cv%  note
----------------------------------------------------------------------------------
```

The CSV contains one row per measured production run. It includes the run ID,
the performance value, aggregate statistics, atom count, tuned kernel values,
input options, and paths to the raw GENESIS logs.

## Input Data

The repository stores input data as compressed archives:

```text
data/<system>.tgz
```

When a selected system is missing from `data/<system>/`, the benchmark extracts
the archive automatically while holding the benchmark lock.

The extracted directories and benchmark outputs are generated files and are not
committed:

```text
data/<system>/
results/
```

## Regenerate Inputs

The checked-in input files live in `inputs/`. Regenerate them with:

```bash
python generate_inputs.py
```

The generated input matrix covers 54 benchmark cells:

- systems: `dhfr, apoa1, uun, factorix, bpti, dppc, ake, stmv, cellulose`
- ensembles: `nve, nvt, npt`
- time steps: `2fs, 4fs` for every system

## Notes

- Benchmark runs are serialized with an advisory lock at `/tmp/bench.lock`.
- The lock file may remain on disk; only a live process holding the lock blocks
  another benchmark.
- The default autotuner is `kernel`, which is the stable option for the current
  benchmark driver.
- `dhfr` is the standard AMBER JAC PME benchmark system with 23,558 atoms. The
  current AMBER suite publishes it as `PME/Topologies/JAC.prmtop` and
  `PME/Coordinates/JAC.inpcrd` in the
  [AMBER20 benchmark suite](https://ambermd.org/Amber20_Benchmark_Suite.tar.gz).
  `data/dhfr.tgz` retains the normal-mass, GENESIS-compatible conversion used
  for the 2 fs input; the 4 fs input applies GENESIS runtime HMR.
- Every rigid system explicitly selects `cons_scheme = MSHAKE` with
  `iter_solute = 3` and `iter_water = 3`.
- Every topology uses normal 1.008 amu hydrogen masses. The 2 fs inputs retain
  those masses, while every 4 fs input enables GENESIS runtime HMR with
  `hmr_target = all`, `hmr_ratio = 3.0`, and
  `hydrogen_mass_upper_bound = 3.3`.
- FactorIX does not load its legacy MD restart because that file contains
  velocities generated for the former pre-HMR topology. GENESIS instead
  initializes velocities after selecting the active 2 fs or 4 fs masses.
- `data/stmv.tgz` is compressed so no single tracked file is larger than
  GitHub's 100 MB file limit.
