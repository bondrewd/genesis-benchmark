# Multi-engine GPU MD benchmark inputs

This repository generates comparable GPU molecular-dynamics inputs for
GENESIS, OpenMM, GROMACS, Amber, and NAMD. Each benchmark system keeps the
single force field and water model supplied by its source archive. The
generator does not reparameterize a system merely to make another engine
available; unsupported engine/model combinations fail closed.

## Physical models

| System | Native model | Water |
| --- | --- | --- |
| DHFR | archive AMBER JAC | TIP3P |
| ApoA1 | CHARMM27 protein/lipid | modified TIP3P |
| UUN | CHARMM36 protein/lipid | modified TIP3P |
| FactorIX | archive AMBER FactorIX model | TIP3P |
| BPTI | Amber03 | TIP3P |
| DPPC | CHARMM36 lipid | modified TIP3P |
| AKE | archive AMBER model | TIP3P |
| STMV | archive AMBER protein/RNA model | TIP3P |
| Cellulose | archive AMBER carbohydrate model | TIP3P |

Every generated explicit-PME system is formally neutral. DHFR and ApoA1 are
neutralized reproducibly by replacing native water molecules with native-family
counterions. Already-neutral source systems are not compositionally changed.

## Protocol

GENESIS inputs use the requested production controls:

- VVER integration;
- BUSSI thermostat and LANGEVIN barostat;
- `thermostat_period`, `barostat_period`, and `baroscale_period` of 10 steps at
  2 fs and 5 steps at 4 fs;
- MSHAKE with `iter_solute = 2`, plus SETTLE for rigid water;
- solute-only HMR at 4 fs, leaving rigid-water masses unchanged;
- periodic center-of-mass motion removal disabled consistently across engines;
- scalar reporting only during timed dynamics, with periodic trajectories and
  restart checkpoints disabled;
- 9/11 A cutoff/pair-list distance for AMBER-family models;
- 12/14 A cutoff/pair-list distance for CHARMM-family models.

Other engines use the closest complete native implementation and record any
non-equivalence directly in the generated input. In particular, GROMACS uses
`md-vv`, v-rescale, C-rescale and LINCS; OpenMM uses CUDA
LangevinMiddle/MonteCarloBarostat; Amber26 uses GPU `pmemd.cuda` with Bussi,
SHAKE and its Monte Carlo barostat; NAMD uses stochastic velocity rescaling,
Langevin piston, SHAKE and SETTLE. No CPU execution path is supported for the
cross-engine comparison. Amber's Bussi implementation has a coupling-time
control but no independent every-N thermostat application interval, so the
requested 10/5-step cadence is recorded as provenance rather than mapped to a
nonexistent Amber option. With `tau_t=5 ps`, Amber applies its coupled Bussi
update every 2,500 steps at 2 fs and every 1,250 steps at 4 fs.

All engines use a PME direct-space tolerance of `1e-5`. GENESIS retains the
established per-system FFT grids; engines whose native setup selects grids from
a target spacing retain that native grid policy.

## Installed engines

The local validated installations are expected at:

```text
GENESIS  ../genesis-mkl-private-gpu/src/spdyn_singlempi/spdyn
OpenMM   /home/diego/miniforge3/envs/OpenMM
GROMACS  /home/diego/gpu-development/gromacs-v2026.3-install/bin/gmx
Amber    /home/diego/miniforge3/envs/Amber26-GPU/bin/pmemd.cuda
NAMD     /home/diego/gpu-development/NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3
```

The Amber environment is intentionally separate from the existing
`AmberTools26` environment. CPU `sander` and QUICK QM/MM CUDA are not accepted
as substitutes for classical GPU `pmemd.cuda`. Its wrapper selects the
Release-26 SPFP binary installed under
`/home/diego/gpu-development/amber/amber26-gpu`; that binary was built for
CUDA 13.1 / `sm_120` and has SHA-256
`31de7efa8680f5b138aced02d075a81efd8325783aeace63117af85c1cd4cf66`.
GROMACS is the clean official `v2026.3` tag at commit
`121090014570a53a17ea391bcddae45e5ea05eb4`.

The installed NAMD 3.0.3 CUDA build is deliberately unavailable for STMV.
Its CUDA tile-list kernel fails for this system, and Compute Sanitizer confirmed
out-of-bounds global reads. The same native STMV model remains available in
GENESIS, OpenMM, GROMACS, and Amber.
Cellulose is supported by NAMD; its large-system launch uses the documented
`+p8 +setcpuaffinity` route on one GPU.

## Prepare native assets

Source archives live in `data/<system>.tgz`. Asset preparation requires OpenMM
and ParmEd but adds no runtime dependency to the generated inputs:

```bash
uv pip install \
  --python /home/diego/miniforge3/envs/OpenMM/bin/python \
  parmed==4.3.1

conda run -n OpenMM python prepare_variants.py --list
conda run -n OpenMM python prepare_variants.py --force
```

Prepared files are reproducible build products under ignored `data/variants/`
directories. Their manifests pin source archives, input members, physical
charges, HMR mass transfers, generated files, and validation results.

## Generate inputs

Show the current native-model support matrix:

```bash
python3 generate_inputs.py --list
```

Generate every supported input:

```bash
python3 generate_inputs.py
```

Inputs are separated by program:

```text
inputs/GENESIS/
inputs/OPENMM/
inputs/GROMACS/
inputs/AMBER/
inputs/NAMD/
```

Filenames are `system__ensemble__timestep.ext`; there is no force-field variant
selector because each system has exactly one native model. Unsuffixed files are
the closest GENESIS-matched execution profile. A distinct engine-native
candidate adds `__native` before the extension.

Only nonredundant profiles are emitted:

| Engine | Generated execution profiles |
| --- | --- |
| GENESIS | matched reference |
| OpenMM | native CUDA route; no built-in Bussi/Langevin-barostat match |
| GROMACS | matched `md-vv`, plus native `md` candidate |
| Amber | one matched route, which is also its native GPU route |
| NAMD | matched route; distinct native Monte Carlo-pressure candidate for NPT |

The `native` label denotes an engine-native/performance candidate, not a
pre-measured winner. Performance conclusions require benchmarking both routes.

Example customization:

```bash
python3 generate_inputs.py \
  --systems dhfr,bpti \
  --engines GENESIS,OPENMM,GROMACS,AMBER,NAMD \
  --ensembles npt \
  --dt 2,4 \
  --nsteps 10000 \
  --output-period 1000 \
  --temperature 300 \
  --pressure 1 \
  --seed 314159 \
  --pair-list-skin 2
```

Explicit engine selections fail if any requested engine cell lacks a fully
validated native representation. Profile-only filters emit the profiles that
exist across engines and report unsupported cells; a default unfiltered
generation likewise skips unsupported cells and reports their count.

## GPU execution

Generated inputs contain the concrete topology and coordinate paths. These
DHFR 4 fs examples show the required GPU launch forms; use the corresponding
paths recorded in another system's input when changing systems.

```bash
flock -x /tmp/bench.lock \
  mpirun -np 1 ../genesis-mkl-private-gpu/src/spdyn_singlempi/spdyn \
  inputs/GENESIS/dhfr__npt__4fs.inp

flock -x /tmp/bench.lock \
  conda run -n OpenMM env GENESIS_BENCHMARK_ROOT="$PWD" \
  python inputs/OPENMM/dhfr__npt__4fs__native.py

flock -x /tmp/bench.lock bash -c '
  gmx=/home/diego/gpu-development/gromacs-v2026.3-install/bin/gmx
  "$gmx" grompp -maxwarn 1 \
    -f inputs/GROMACS/dhfr__npt__4fs.mdp \
    -c data/variants/dhfr/native_amber_jac_tip3p/system.g96 \
    -p data/variants/dhfr/native_amber_jac_tip3p/system_hmr.top \
    -po /tmp/dhfr-mdout.mdp \
    -o /tmp/dhfr.tpr
  "$gmx" mdrun -s /tmp/dhfr.tpr -deffnm /tmp/dhfr \
    -ntmpi 1 -ntomp 1 -nb gpu -pme gpu -pmefft gpu -bonded gpu -update cpu \
    -tunepme no -gpu_id 0 -noconfout -cpt -1
'

flock -x /tmp/bench.lock bash -c '
  gmx=/home/diego/gpu-development/gromacs-v2026.3-install/bin/gmx
  "$gmx" grompp -maxwarn 1 \
    -f inputs/GROMACS/dhfr__npt__4fs__native.mdp \
    -c data/variants/dhfr/native_amber_jac_tip3p/system.g96 \
    -p data/variants/dhfr/native_amber_jac_tip3p/system_hmr.top \
    -po /tmp/dhfr-native-mdout.mdp \
    -o /tmp/dhfr-native.tpr
  "$gmx" mdrun -s /tmp/dhfr-native.tpr -deffnm /tmp/dhfr-native \
    -ntmpi 1 -ntomp 1 -nb gpu -pme gpu -pmefft gpu -bonded gpu \
    -update gpu -gpu_id 0 -noconfout -cpt -1
'

flock -x /tmp/bench.lock \
  conda run -n Amber26-GPU pmemd.cuda -O \
  -i inputs/AMBER/dhfr__npt__4fs.in \
  -p data/variants/dhfr/native_amber_jac_tip3p/system_hmr.prmtop \
  -c data/variants/dhfr/native_amber_jac_tip3p/system.inpcrd \
  -o /tmp/dhfr.mdout -r /tmp/dhfr.restrt

flock -x /tmp/bench.lock \
  /home/diego/gpu-development/NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3 \
  +p1 +devices 0 inputs/NAMD/dhfr__npt__4fs.namd

# NAMD's distinct native NPT candidate uses Monte Carlo pressure control.
flock -x /tmp/bench.lock \
  /home/diego/gpu-development/NAMD_3.0.3_Linux-x86_64-multicore-CUDA/namd3 \
  +p1 +devices 0 inputs/NAMD/dhfr__npt__4fs__native.namd
```

GROMACS 2026.3 does not support GPU update/constraint offload with `md-vv`, so
the matched profile keeps its update stage on the CPU while nonbonded, PME, and
supported bonded work are explicitly offloaded. The separately labeled native
profile uses `md`, engine-managed coupling/buffering, and GPU update so the
performance/algorithm tradeoff can be measured rather than assumed. For native
NVE, the generated Verlet-buffer tolerance is capped by temperature and run
duration so GROMACS's estimated accumulated drift remains at or below 1%.
The single allowed `grompp` warning is the expected notice that periodic
center-of-mass removal is deliberately disabled; automation rejects any other
warning or a different warning count.

The GENESIS benchmark driver additionally rejects logs that do not prove the
real-space, reciprocal-space, pair-list, and nonbonded GPU routes.

## GENESIS measurements

`run_benchmark.py` remains the focused GENESIS GPU measurement driver. It reads
only `inputs/GENESIS`, verifies immutable scientific controls against a fresh
canonical render, serializes GPU access through `/tmp/bench.lock`, and retains
raw logs plus CSV results.

```bash
python3 run_benchmark.py \
  --spdyn ../genesis-mkl-private-gpu/src/spdyn_singlempi/spdyn \
  --systems dhfr \
  --ensembles npt \
  --dt 4 \
  --warmup 1 \
  --measure 5
```

Short acceptance runs are not performance measurements: one 10,000-step GPU
run is used per advertised engine/system/execution-profile combination, reduced
to 1,000 steps for Cellulose and STMV. Repeated warmups and retained samples are
required only when reporting comparative performance.
