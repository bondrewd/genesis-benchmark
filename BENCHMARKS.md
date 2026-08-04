# Benchmarks

Historical entries that use the system name `dhfr` refer to the current
`dhfr_27k` dataset.

## 2026-07-08 19:18 JST - DHFR 4 fs benchmark-driver fix

Repository commit: `b70d450` (`genesis-benchmark`, dirty with this change).
Benchmark patch diff SHA256 before this append-only log update:
`30d96199d527971679a4149782c28adb5df9cfb4e8de025b0ab09dd6aae8689c`
for `README.md`, `generate_inputs.py`, `run_benchmark.py`, and `inputs/*.inp`.
GENESIS commit: `dbea5e7f0` (`/home/diego/Repos/genesis-mkl-private`, dirty before this task).
GENESIS dirty diff SHA256:
`c9f17da64c872826fa3ba903478dcdfb873bba3bea67c75972d08cfb4ed49994`.
GENESIS dirty files:
`src/gpu/gpu_boundary.cu`, `src/spdyn_singlempi/sp_boundary.fpp`,
`src/spdyn_singlempi/sp_enefunc.fpp`, `src/spdyn_singlempi/sp_enefunc_str.fpp`,
`src/spdyn_singlempi/sp_energy.fpp`, `src/spdyn_singlempi/sp_md_vverlet.fpp`,
`src/spdyn_singlempi/sp_setup_spdyn.fpp`, `src/spdyn_singlempi/sp_update_domain.fpp`.
GENESIS dirty diffstat: 8 files changed, 417 insertions, 332 deletions.
Binary: `/home/diego/Repos/genesis-mkl-private/src/spdyn_singlempi/spdyn`.
Environment: `GENESIS_GPU_PROFILE=0`, `OMP_NUM_THREADS=1`, `HWLOC_COMPONENTS=x86`, advisory lock `/tmp/bench.lock`.
Hardware: Intel Core Ultra 9 275HX, NVIDIA GeForce RTX 5070 Ti Laptop GPU, driver 595.71.05.
Python: 3.12.13 via `uv run`.

### Pre-fix reproducer: DHFR/NVT 1000 steps

Command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nvt --dt 2,4 --warmup 1 --measure 5 --nsteps 1000 --tune-nsteps 1000 --timeout 300 --timestamp dhfr-repro-nvt-current
```

Result: 4 fs was unstable, with 102.21% CV and alternating slow/fast timings.

```text
dhfr,nvt,2fs,573.874,0.412,0.07,5,raw=573.350|573.950|574.320|573.560|574.190
dhfr,nvt,4fs,779.164,796.391,102.21,5,raw=197.530|201.440|1645.770|1657.330|193.750
```

### Post-fix acceptance: DHFR all ensembles 10000 steps

Command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nve,nvt,npt --dt 2,4 --warmup 2 --measure 5 --nsteps 10000 --tune-nsteps 10000 --timeout 600 --timestamp dhfr-all-final-10k
```

Result: all selected cells completed with no errors. The 4 fs ns/day is approximately 2x the matching 2 fs cell.
Pinned inputs from `results/work`:

```text
dhfr_nve_2fs.pinned.inp sha256=27404c32dfee48118be3e881623c6059509b3950d5b59143a475629b00de8dbc kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
dhfr_nve_4fs.pinned.inp sha256=a062fc94b8dbe29fbc49b9290808c45679cd0bde2f2c1b14737b2dafcabb0469 kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
dhfr_nvt_2fs.pinned.inp sha256=74ce206b1295fea7139f1c96c48fc8ee8c41bc7b33a146c0166d06c814031a5b kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
dhfr_nvt_4fs.pinned.inp sha256=15c0eebee58e174682ab42a7646f859a59d0324f18505154652ad47b8c30ed55 kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
dhfr_npt_2fs.pinned.inp sha256=20cbc8303715aa19533589c89eecf9766a780054d9ad34af7d73bffdd3fd8b7d kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
dhfr_npt_4fs.pinned.inp sha256=bb1ec84fa4e83db769c82ad1f3e3d17542b571da4da551810889c5308d92bf71 kernel=(spread=128,influence=256,bonded=128,constraints=128,nonbond_inter_mb=8,nonbond_intra_mb=6)
```

```text
dhfr,nve,2fs,701.028,1.849,0.26,5,raw=703.100|701.970|698.140|700.690|701.240
dhfr,nve,4fs,1397.328,1.103,0.08,5,raw=1398.530|1398.350|1397.090|1395.920|1396.750,ratio=1.993
dhfr,nvt,2fs,671.782,0.143,0.02,5,raw=671.630|671.710|671.990|671.860|671.720
dhfr,nvt,4fs,1339.402,1.768,0.13,5,raw=1337.580|1340.970|1337.860|1341.470|1339.130,ratio=1.994
dhfr,npt,2fs,579.310,0.652,0.11,5,raw=579.450|579.870|579.830|579.120|578.280
dhfr,npt,4fs,1174.072,1.074,0.09,5,raw=1173.830|1174.640|1173.660|1172.690|1175.540,ratio=2.027
```

Coverage note: this timed acceptance run covers DHFR only. The same runtime-HMR
input-generation policy is applied to other non-pre-HMR 4 fs systems in this
change, but their timed 2 fs/4 fs ratios were not measured in this entry. Run the
full selected matrix before treating the whole suite as release-validated.

## 2026-07-08 20:12 JST - data archive packaging checks

Repository commit: `b70d450` (`genesis-benchmark`, dirty with this change).
Scope: packaging/non-GPU verification for storing committed inputs as
`data/*.tgz` and extracting `data/<system>/` on demand. No performance timing
claim is made by this entry.

Commands and results:

```bash
for f in data/*.tgz; do tar -tzf "$f" >/dev/null || exit 1; done
# passed for all eight archives

python3 -c 'import run_benchmark; print(run_benchmark.ensure_system_data(["dhfr"]))'
# ['dhfr'] on a missing data/dhfr directory

python3 -c 'import run_benchmark; print(run_benchmark.ensure_system_data(["dhfr"]))'
# [] after a complete extraction with matching .archive_sha256

rm data/dhfr/equil.rst
python3 -c 'import run_benchmark; print(run_benchmark.ensure_system_data(["dhfr"]))'
# ['dhfr']; data/dhfr/equil.rst restored at 1338827 bytes

find . -type f -size +100M -printf '%p %s\n'
# no output

git ls-files -s | awk '{print $4}' | while read -r f; do
  if [ -f "$f" ]; then
    s=$(stat -c%s "$f")
    if [ "$s" -gt 100000000 ]; then printf '%s %s\n' "$s" "$f"; fi
  fi
done
# no output
```

## 2026-07-09 09:44 JST - pre-HMR 2 fs setup regression check

Repository commit: `5ae4972` (`genesis-benchmark`, dirty with this change).
Scope: verify that FactorIX and STMV 2 fs inputs with pre-HMR topologies pass
GENESIS setup and run a tiny 10-step window after adding
`hydrogen_mass_upper_bound = 3.3`. This is not a performance baseline.

Command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems factorix,stmv --ensembles nve,nvt,npt --dt 2 --tune kernel --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --timeout 300 --timestamp verify-prehmr-2fs
```

Result: all six cells completed with no setup/tune/measure failures.

```text
factorix,nve,2fs,86.27,0.00,0.0,1,raw=86.270
factorix,nvt,2fs,85.46,0.00,0.0,1,raw=85.460
factorix,npt,2fs,86.35,0.00,0.0,1,raw=86.350
stmv,nve,2fs,2.58,0.00,0.0,1,raw=2.580
stmv,nvt,2fs,2.59,0.00,0.0,1,raw=2.590
stmv,npt,2fs,2.47,0.00,0.0,1,raw=2.470
```

## 2026-07-09 11:17 JST - per-run logging and CSV smoke check

Repository commit: `5ae4972` (`genesis-benchmark`, dirty with this change).
Scope: functional check for the benchmark driver's per-result log bundle and
per-measurement CSV rows. This is a 10-step smoke test, not a performance
baseline.

Command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nve --dt 2 --tune kernel --warmup 1 --measure 2 --nsteps 10 --tune-nsteps 10 --timeout 120 --timestamp log-smoke2
```

Result: completed without failures. The driver wrote
`results/log-smoke2.csv`, `results/log-smoke2/summary.log`, exact generated
inputs under `results/log-smoke2/inputs/`, one autotune log, one warmup log, and
two measured production logs. The CSV contains two rows with `run_id`,
`ns_per_day`, `num_atoms=27346`, tuned kernel columns, input option columns, and
log-path columns.

```text
dhfr,nve,2fs,mean=208.795,median=208.795,std=0.007,cv=0.00,n=2,raw=208.790|208.800
```

## 2026-07-09 12:01 JST - benchmark.log / summary.log split smoke check

Repository commit: `5ae4972` (`genesis-benchmark`, dirty with this change).
Scope: functional check that the full progress transcript is written to
`benchmark.log` and `summary.log` contains only the final aggregate table. This
is a 10-step smoke test, not a performance baseline.

Command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nve --dt 2 --tune kernel --warmup 1 --measure 2 --nsteps 10 --tune-nsteps 10 --timeout 120 --timestamp log-summary-smoke
```

Result: completed without failures. `results/log-summary-smoke/benchmark.log`
contains the full benchmark progress log. `results/log-summary-smoke/summary.log`
starts with the `=== ns/day ... ===` table header and excludes setup/progress
lines and the `CSV written` footer.

```text
dhfr,nve,2fs,mean=210.605,median=210.605,std=0.870,cv=0.41,n=2,raw=211.220|209.990
```

## 2026-07-09 12:13 JST - eneout_period CLI smoke check

Repository commit: `5ae4972` (`genesis-benchmark`, dirty with this change).
Scope: functional check that the driver accepts `--eneout-period`, validates the
step-count divisibility rule, and writes the requested value into the generated
tune and pinned inputs. This is a 10-step smoke test, not a performance
baseline.

Validation checks:

```bash
python3 run_benchmark.py --nsteps 1500 --tune-nsteps 50000 --eneout-period 1000
# exits 2: --nsteps (1500) must be a multiple of --eneout-period (1000)

python3 run_benchmark.py --nsteps 100000 --tune-nsteps 50500 --eneout-period 1000
# exits 2: --tune-nsteps (50500) must be a multiple of --eneout-period (1000)
```

Smoke command:

```bash
uv run run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nve --dt 2 --tune kernel --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 120 --timestamp eneout-smoke
```

Result: completed without failures. Both
`results/eneout-smoke/inputs/dhfr_nve_2fs.tune.inp` and
`results/eneout-smoke/inputs/dhfr_nve_2fs.pinned.inp` contain
`nsteps = 10` and `eneout_period = 10`; the CSV row contains
`input_dynamics_nsteps=10` and `input_dynamics_eneout_period=10`.

```text
dhfr,nve,2fs,mean=201.820,median=201.820,std=0.000,cv=0.00,n=1,raw=201.820
```

## 2026-07-09 13:09 JST - code audit smoke check

Repository commit: `8a72c94` (`genesis-benchmark`, dirty with this change).
Scope: functional check after renaming Python variables and simplifying comments.
This is a 10-step smoke test, not a performance baseline.

Command:

```bash
python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr --ensembles nve --dt 2 --tune kernel --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 120 --timestamp audit-smoke
```

Result: completed without failures. The CSV contains one measured row with
`input_dynamics_nsteps=10` and `input_dynamics_eneout_period=10`; `summary.log`
contains only the aggregate table and `benchmark.log` contains the full progress
log.

```text
dhfr,nve,2fs,mean=205.500,median=205.500,std=0.000,cv=0.00,n=1,raw=205.500
```

## 2026-07-09 14:10 JST - Amber PME import smoke checks

Repository commit: `2940e00` (`genesis-benchmark`, dirty with this change).
Host: `legion`; GPU: NVIDIA GeForce RTX 5070 Ti Laptop GPU, driver 595.71.05.
Scope: functional smoke checks for the Amber PME additions and the DHFR rename.
These are 10-step checks with no autotune, not performance baselines.

Commands:

```bash
python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr_23k --ensembles nve,nvt,npt --dt 2,4 --tune none --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 300 --timestamp pme-dhfr23-runtime-hmr-smoke

python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems cellulose --ensembles nve,npt --dt 2,4 --tune none --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 600 --timestamp pme-cellulose-smoke

python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems cellulose --ensembles npt --dt 2,4 --tune none --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 600 --timestamp pme-cellulose-npt-smoke2

python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems cellulose --ensembles nvt --dt 2,4 --tune none --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 600 --timestamp pme-cellulose-nvt-smoke

python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn --systems dhfr_27k --ensembles nve --dt 2 --tune none --warmup 0 --measure 1 --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 300 --timestamp pme-dhfr27-rename-smoke
```

Results: all final smoke checks completed without failures. Amber's old-format
JAC 2 fs topology failed before conversion to modern Amber topology format, so
`data/dhfr_23k/prmtop` was converted with ParmEd while preserving 1.008 amu
hydrogen masses. DHFR 23k uses that normal-mass topology for both 2 fs and 4 fs;
the 4 fs inputs use GENESIS runtime HMR.
The first cellulose smoke run passed NVE but failed NPT without explicit
`box_size_*` values; `pme-cellulose-npt-smoke2` supersedes those failed NPT
cells after adding the cellulose box dimensions to the generated NPT inputs.

```text
dhfr_23k,nve,2fs,atoms=23558,raw=98.270
dhfr_23k,nve,4fs,atoms=23558,raw=187.850
dhfr_23k,nvt,2fs,atoms=23558,raw=98.590
dhfr_23k,nvt,4fs,atoms=23558,raw=195.290
dhfr_23k,npt,2fs,atoms=23558,raw=97.590
dhfr_23k,npt,4fs,atoms=23558,raw=189.050
cellulose,nve,2fs,atoms=408609,raw=7.080
cellulose,nve,4fs,atoms=408609,raw=14.040
cellulose,npt,2fs,atoms=408609,raw=6.790
cellulose,npt,4fs,atoms=408609,raw=13.330
cellulose,nvt,2fs,atoms=408609,raw=6.750
cellulose,nvt,4fs,atoms=408609,raw=13.960
dhfr_27k,nve,2fs,atoms=27346,raw=206.680
```

## 2026-08-04 - input compatibility and canonical DHFR refresh

Repository base commit: `95152844deacaaaa83360e02bece9d3bc99e4d82`
(`genesis-benchmark`, dirty with this change). GENESIS commit:
`f8de4a3cc3a85d5b56cf852e9bf7850f41872522` (clean).

Scope: input/data correctness only. No performance comparison is claimed.
The active DHFR system is now the standard 23,558-atom JAC PME benchmark, and
the 27,346-atom dataset is removed from the active matrix. The AMBER20 suite
downloaded from `https://ambermd.org/Amber20_Benchmark_Suite.tar.gz` has SHA256
`7fee3a02f85f0eb1d07b6eb8902516947f29774ac870dbe17a7cb134d916f401` and
contains `PME/Topologies/JAC.prmtop` with `NATOM = 23558`. The retained
normal-mass GENESIS conversion has the same atom names, atom types, charges,
bond count, angle count, dihedral count, residue count, atom count, and total
mass as that current JAC topology; its normal hydrogen masses permit a genuine
normal-mass 2 fs input and GENESIS runtime HMR at 4 fs.

The packaged FactorIX and STMV AMBER topologies were also normalized from
pre-applied HMR. The transformation changed only the `%FLAG MASS` section:
each 3.024 amu hydrogen was restored to 1.008 amu and 2.016 amu was returned
to its single bonded heavy atom. FactorIX restored 2,817 hydrogens bonded to
1,851 heavy atoms; STMV restored 77,849 hydrogens bonded to 52,761 heavy atoms.
Total topology mass was preserved at 554,212.334049 amu and 6,695,311.430000
amu, respectively. All packaged systems now use normal hydrogen masses, with
GENESIS runtime HMR enabled only by the 4 fs inputs.

The FactorIX MD restart is retained in the source archive for provenance but
is no longer referenced by benchmark inputs. It contains velocities generated
for the former pre-HMR topology. Omitting it makes GENESIS initialize fresh
velocities after runtime mass selection, avoiding an HMR/normal-mass kinetic
energy mismatch in the new 2 fs cells.

Canonical local archive hashes:

```text
data/dhfr.tgz    15491921bb5ccfae41d3f8ea37270fd8206180a82aca17a760a5d952ac896b6f
dhfr/prmtop      027f7be194c2ad285281528011341ec40ec050f7b3acdf487094d7e2284993bf
dhfr/inpcrd      7f084ccce438dc61e49bbf6c3ee2d2a56cf1d057dcd64cace42ae7b2a0105cfe
data/factorix.tgz d2d7950bac0c384fbc84ccc00845e831e7b73750640ce78d0fec2b5f73d653cc
factorix/FactorIX.prmtop 14c4527b64c5f4b6bef1fbbd90fec0fa6b8282bcc5517d2c0efae676e74d3838
data/stmv.tgz    d991da1c1ee17f442730c41a68c9d3adc54b355ba1dc2d7e748b7f826bc7a735
stmv/prmtop      94ccc723db70b7482f293632f880553ccb3ea73347c903cafd3cc7aecb700068
```

Static verification:

```bash
python3 -m unittest -v
# 9 tests passed: exact 54-cell matrix, retired-key exclusion, global equals
# alignment, MSHAKE/3/3 constraints, runtime-HMR policy, packaged AMBER/PSF/
# GROMACS mass normalization, current kernel tuner parsing, and DHFR archive
# NATOM/prefix checks.

python3 -m py_compile generate_inputs.py run_benchmark.py test_benchmark_inputs.py
# passed

cd ../genesis-mkl-private && \
  python3 tests/unit_tests/test_constraint_scheme_mapping.py
# constraint scheme mapping: 9 passed, 0 failed

git diff --check
# passed
```

A lock-serialized all-54-cell 10-step GENESIS smoke was attempted with:

```bash
python3 run_benchmark.py --spdyn ../genesis-mkl-private/src/spdyn_singlempi/spdyn \
  --ensembles nve,nvt,npt --dt 2,4 --tune none --warmup 0 --measure 1 \
  --nsteps 10 --tune-nsteps 10 --eneout-period 10 --timeout 600 \
  --timestamp input-refresh-smoke-native
```

The current host exposed no CUDA device (`nvidia-smi` could not communicate
with the driver), so GENESIS stopped at `gpu_info.cu` before input parsing with
`no CUDA-capable device is detected`. Consequently, no integration or
performance result from that attempted run is counted as a pass.
