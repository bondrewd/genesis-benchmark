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
  --systems dhfr_27k \
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

- all systems: `dhfr_27k,dhfr_23k,apoa1,uun,factorix,bpti,dppc,ake,stmv,cellulose`
- all ensembles: `nve,nvt,npt`
- both time steps: `2fs,4fs`, except `factorix` and `stmv`, which run only `4fs`
- kernel autotuning
- `100000` production steps
- `50000` autotuning steps
- `eneout_period = 1000`
- `1` MPI process
- `1` OpenMP thread

## Common Examples

Run only the original 27k-atom DHFR NVE at both 2 fs and 4 fs:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr_27k \
  --ensembles nve \
  --dt 2,4
```

Run Amber's 23k-atom DHFR and cellulose systems:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr_23k,cellulose \
  --ensembles nve \
  --dt 2
```

Use shorter runs for a smoke test:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr_27k \
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
  --systems dhfr_27k \
  --mpi-procs 1 \
  --omp-threads 1
```

Override the energy-output period:

```bash
python run_benchmark.py \
  --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --systems dhfr_27k \
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
dhfr_27k   nve   2fs        662.52       662.47       0.89    0.1%
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

- systems: `dhfr_27k, dhfr_23k, apoa1, uun, factorix, bpti, dppc, ake, stmv, cellulose`
- ensembles: `nve, nvt, npt`
- time steps: `2fs, 4fs` for normal-mass topologies, and only `4fs` for
  `factorix` and `stmv`

## Notes

- Benchmark runs are serialized with an advisory lock at `/tmp/bench.lock`.
- The lock file may remain on disk; only a live process holding the lock blocks
  another benchmark.
- The default autotuner is `kernel`, which is the stable option for the current
  benchmark driver.
- The original DHFR system is named `dhfr_27k`; Amber's PME DHFR system is
  named `dhfr_23k`.
- Normal 1.008 amu hydrogen topologies are used for 2 fs and can use GENESIS
  runtime HMR for 4 fs. Topologies that already contain 3.024 amu hydrogens are
  treated as HMR-only inputs; those inputs set `hydrogen_mass_upper_bound = 3.3`
  and do not use the `hydrogen_mr`, `hmr_target`, or `hmr_ratio` runtime-scaling
  options.
- `factorix` and `stmv` use HMR topologies, so the benchmark does not generate
  or run 2 fs inputs for those systems.
- `data/stmv.tgz` is compressed so no single tracked file is larger than
  GitHub's 100 MB file limit.
