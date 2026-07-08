# Benchmarks

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
