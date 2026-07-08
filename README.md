# GENESIS spdyn GPU performance benchmark

Self-contained GPU benchmarking harness for GENESIS `spdyn`. It measures the
profiler-independent **`[PERFORMANCE]` ns/day** for a matrix of

**8 systems x {NVE, NVT, NPT} x {2 fs, 4 fs} = 48 input files**,

each first *tuned* (autotuners on), then *pinned* (tuned values hard-coded,
autotuners off), then *warmed up* and *measured*.

Everything the runs need is inside this directory. Git stores one compressed
archive per system (`data/<sys>.tgz`); `run_benchmark.py` extracts
`data/<sys>/` on demand when a selected system directory is missing. The suite
does not depend on `tests/` or any loose directory in the repo.

---

## Directory layout

```
benchmark/
  README.md               # this file
  generate_inputs.py      # regenerates the 48 inputs/ files from per-system templates
  run_benchmark.py        # the driver: tune -> pin -> warm-up -> measure
  data/<sys>.tgz          # committed topology + coordinate archive for each system
  data/<sys>/             # extracted local topology + coordinates (generated, ignored)
  inputs/<sys>_<ens>_<dt>.inp   # the 48 base inputs (referenced by the driver)
  results/<timestamp>.csv       # one local CSV per driver run (generated, ignored)
  results/work/                 # generated .tune.inp / .pinned.inp per cell (ignored)
```

### The 8 systems (`<sys>`)

| `<sys>`    | benchmark | force field | source (tests/performance_tests reference) |
|------------|-----------|-------------|---------------------------------------------|
| `dhfr`     | DHFR / JAC        | AMBER    | 04/05/06 (`jac_amber`)  |
| `apoa1`    | ApoA1             | CHARMM   | 01/02/03 (`apoa1`)      |
| `uun`      | UUN               | CHARMM   | 07/08/09 (`uun`)        |
| `factorix` | Factor IX         | AMBER    | 25 (NPT; `factorix`)    |
| `bpti`     | BPTI              | GROAMBER | 16/17/18 (`bpti`)       |
| `dppc`     | DPPC bilayer      | CHARMM   | 10/11/12 (`dppc`)       |
| `ake`      | adenylate kinase  | AMBER    | 13/14/15 (`ake`)        |
| `stmv`     | STMV (~1M atoms)  | AMBER    | loose `STMV_production_NPT_4fs/` |

Force-field / PME / cutoff / box parameters for every cell come from the matching
reference input above. `factorix` and `stmv` only had NPT references; their NVE and
NVT variants are derived by dropping the barostat (NPT->NVT) and the thermostat
(NVT->NVE).

### The 48-cell matrix (`<ens>` x `<dt>`)

* **Ensemble** (`nve` / `nvt` / `npt`), uniform across the matrix so cells are comparable:
  * `nve` : `tpcontrol = NO`
  * `nvt` : `tpcontrol = BUSSI`, `thermostat_period = 10`
  * `npt` : `tpcontrol = BUSSI`, `thermostat_period = 10`, `barostat_period = 10`, `pressure = 1.0`
  * all use `group_tp = YES`, `temperature = 300`, fixed `iseed = 314159` for reproducibility.
* **Timestep** (`2fs` / `4fs`):
  * `2fs` : `timestep = 0.002`, `rigid_bond = YES`
  * `4fs` : `timestep = 0.004`, `rigid_bond = YES` **+ hydrogen-mass repartitioning**.
    Inputs whose topology already contains redistributed HMR masses (Factor IX,
    STMV) only raise `hydrogen_mass_upper_bound = 3.3` so GENESIS still
    recognizes those heavier H atoms as hydrogens. Inputs with ordinary hydrogen
    masses enable GENESIS runtime HMR with `hydrogen_mr = YES`, `hmr_target =
    all`, and `hmr_ratio = 3.0`, plus the same 3.3 recognition threshold.

The base inputs ship with `nsteps = 10000`, `eneout_period = 1000` (so they run
stand-alone), but the driver overrides `nsteps`/`eneout_period` for every run.
`[OUTPUT]` is empty, so no trajectory/restart files are written.

### Input data archives

The committed input data are the eight per-system archives under `data/`:

```
data/ake.tgz      data/apoa1.tgz   data/bpti.tgz   data/dhfr.tgz
data/dppc.tgz     data/factorix.tgz data/stmv.tgz  data/uun.tgz
```

On a fresh checkout, the first benchmark run for a selected system extracts only
that system archive while holding the benchmark lock. The uncompressed
`data/<sys>/` directories are generated files and are ignored by Git. You may
remove them locally at any time; the next run recreates the directories from the
archives.

### STMV data (large)

STMV's uncompressed `prmtop` (203 MB) + `inpcrd` (78 MB) are packaged into
`data/stmv.tgz` (about 48 MB) so the repository can be pushed without a single
file over 100 MB. STMV is ~1M atoms; expect it to need a large GPU and to be
very slow. It is **not** part of the quick smoke tests.

If a previous commit already tracked the 203 MB uncompressed `data/stmv/prmtop`,
GitHub can still reject a push because the large blob remains in history. The
working tree fix is to track only `data/*.tgz` and untrack `data/<sys>/`, but an
already-published or already-created history may also need to be rewritten
(`git filter-repo`, an orphan branch, or Git LFS) before pushing to GitHub.

---

## Measurement protocol (per selected `system x ensemble x dt` cell)

`run_benchmark.py` executes, for each cell:

1. **TUNE** — run once with the requested autotuners ON, then parse the
   `[AUTOTUNE]` report to extract the tuned values:
   * kernel block sizes (the `kernel_*` paste-ready lines),
   * `cell_size` (the `Selected configuration: cell_size = X` line),
   * `pairlistdist` + `nbupdate_period` (the neighbor-list candidate marked `(selected)`).
2. **PIN** — write a *fresh* input that hard-codes the tuned values
   (`kernel_*` into `[GPU]`, `cell_size` + `pairlistdist` into `[ENERGY]`,
   `nbupdate_period` into `[DYNAMICS]`) and turns **all** autotuners OFF, so the
   measured runs carry no autotune instrumentation.
3. **WARM-UP** — run the pinned input `--warmup` times (default 2), timings discarded
   (heat the GPU / stabilise clocks).
4. **MEASURE** — run the pinned input `--measure` times (default 10), collect the
   `[PERFORMANCE]` ns/day of each.
5. Report **mean ± std (and cv%)** and append a row to `results/<timestamp>.csv`.

Every `spdyn` launch is
```
GENESIS_GPU_PROFILE=0 OMP_NUM_THREADS=1 HWLOC_COMPONENTS=x86 \
    mpirun -np 1 <repo>/src/spdyn_singlempi/spdyn <input>
```
run from the benchmark root, and the ns/day is read from the `[PERFORMANCE]` line
(GPU-synced, profiler-independent). Runs are **serialised through an advisory
lockfile** (`/tmp/bench.lock`) so concurrent agents never benchmark at the same
time. The lock file itself may remain on disk; this is harmless because the lock
is owned by the live process, not by file existence.

---

## CLI

```
python3 run_benchmark.py [options]

--systems   dhfr,apoa1,...   (default: all 8)
--ensembles nve,nvt,npt      (default: all 3)
--dt        2,4              (default: both, fs)
--warmup    N                (default 2)   discarded warm-up runs
--measure   M                (default 10)  timed runs
--tune      kernel|cell|nblist|all|none  (default: kernel)   comma-list allowed
--full-autotune              shortcut for --tune kernel,cell,nblist
--nsteps        N            (default 10000) measurement-run nsteps (eneout_period is matched)
--tune-nsteps   N            (default 10000) tuning-run nsteps
--timeout    S               (default 7200)  per-run timeout
--lock       PATH            (default /tmp/bench.lock)
--allow-failures             write partial CSV and exit 0 even if cells fail
--timestamp  STR             CSV name stamp (default: local date-time)
--out        PATH            explicit CSV path (overrides --timestamp)
```

### Examples

```bash
# Quick end-to-end check on two fast cells (what was used to validate this harness):
python3 run_benchmark.py --systems dhfr,apoa1 --ensembles nve,nvt --dt 2,4 \
    --tune kernel --warmup 1 --measure 3 --nsteps 2500 --tune-nsteps 3500

# The default kernel-only tune, all 8 systems, full 10-sample measurement:
python3 run_benchmark.py

# Turn on all three autotuners (see the caveat below):
python3 run_benchmark.py --systems ake,dppc --full-autotune

# Neighbor-list tuning only, NVT/NPT only:
python3 run_benchmark.py --tune nblist --ensembles nvt,npt
```

### Regenerating the 48 inputs

```bash
python3 generate_inputs.py     # rewrites inputs/*.inp from the per-system templates
```

---

## Autotuner caveat (important)

The **default is `--tune kernel`**, and kernel-block-size autotuning is robust on
every system tested — it never re-decomposes the cell grid.

The **`cell` and `nblist` autotuners re-decompose the cell grid at runtime**, and in
the current build that path is fragile:

* On systems whose `[BOUNDARY]` specifies a `box_size` (e.g. `dhfr`, `apoa1`) the
  re-decomposition can abort with a GENESIS **cell-overflow**
  (`gpu_domain.cu: cell too large ... exceeds the 256 per-cell capacity`).
* `cell_size_autotune` frequently reports a degenerate `cell_grid = 0 x 0 x 0`
  (the golden-section search does not actually evaluate), i.e. it pins the default.
* On no-`box_size` systems (`ake`, `dppc`) `nblist` tuning works and selects a real
  `pairlistdist`.

This is a GENESIS-side limitation, not a harness bug. The driver handles the
known `cell`/`nblist` re-decomposition overflow gracefully: **if that specific
cell-overflow abort occurs, the cell falls back to a default-pinned measurement**
(all autotuners off, no pinned values) so you still get a valid ns/day, and the
CSV/table `note` column records `tune-failed(cell overflow)`. Other tuning
failures, including kernel-only tuning failures, fail the selected cell. For
dependable end-to-end numbers, use the default `--tune kernel`.

---

## Output

Console prints a table:

```
system     ens   dt         ns/day      +-std     cv%  note
dhfr       nve   2fs        702.45       0.55    0.1%
apoa1      nve   2fs        135.57       0.20    0.1%
```

`results/<timestamp>.csv` columns:
`system, ensemble, dt, ns_per_day_mean, ns_per_day_std, cv_pct, n_measure,
cell_size, pairlistdist, nbupdate_period, tuners, note, raw_ns_per_day`
(the last column is the `|`-joined per-run ns/day so you can inspect variance).

A high `cv%` flags an unstable measurement (thermal throttling or contention on
the shared machine) — re-run that cell when the machine is quiet.
