#!/usr/bin/env python3
"""Generate the 48 benchmark input files: 8 systems x {nve,nvt,npt} x {2fs,4fs}.

Each generated .inp lives in benchmark/inputs/ and references its system's data
with a path RELATIVE TO THE BENCHMARK ROOT (data/<sys>/...). The benchmark driver
(run_benchmark.py) always launches spdyn with cwd = the benchmark root, so those
relative paths resolve correctly.

Base parameters (force field, PME grid, cutoffs, box) are taken from the matching
tests/performance_tests/NN/gpu.inp reference for that system+ensemble. Where a system
only had an NPT reference (factorix, stmv) the NVE/NVT variants are derived by
dropping the barostat (NPT->NVT) and the thermostat (NVT->NVE).

Ensemble policy (uniform across the matrix, so the 48 cells are comparable):
  NVE : tpcontrol = NO
  NVT : tpcontrol = BUSSI, thermostat_period = 10 at 2fs and 5 at 4fs
  NPT : tpcontrol = BUSSI, thermostat_period/barostat_period/baroscale_period
        = 10 at 2fs and 5 at 4fs, pressure = 1.0
Timestep policy:
  2fs : timestep = 0.002, rigid_bond = YES
  4fs : timestep = 0.004, rigid_bond = YES + HMR
        - if the topology already has redistributed HMR masses, only raise
          hydrogen_mass_upper_bound to 3.3 for hydrogen recognition
        - otherwise enable runtime hydrogen_mr with hmr_target = all
  Any topology that already contains redistributed HMR masses uses
  hydrogen_mass_upper_bound = 3.3 at both 2fs and 4fs so GENESIS still
  recognizes those atoms as hydrogens for constraints.
Neighbor-list policy:
  nbupdate_period = 10

Run this from anywhere:  python3 benchmark/generate_inputs.py
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "inputs")

# Per-system definitions. Paths are relative to the benchmark root (data/<sys>/...).
SYSTEMS = {
    # ---- AMBER: DHFR / JAC ----
    "dhfr": dict(
        ff="AMBER",
        input=[
            "prmtopfile = data/dhfr/step3_input.parm7",
            "ambcrdfile = data/dhfr/step3_input.rst7",
            "rstfile    = data/dhfr/equil.rst",
        ],
        energy=[
            "switchdist       = 8",
            "cutoffdist       = 8",
            "pairlistdist     = 9.5",
            "pme_ngrid_x      = 48",
            "pme_ngrid_y      = 48",
            "pme_ngrid_z      = 48",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=["water_model      = WAT"],
        box=("1", "1", "1"), domain=False, npt_drops_box=False,
    ),
    # ---- CHARMM: ApoA1 ----
    "apoa1": dict(
        ff="CHARMM",
        input=[
            "topfile = data/apoa1/top_all27_prot_lipid.rtf",
            "parfile = data/apoa1/par_all27_prot_lipid.prm",
            "psffile = data/apoa1/apoa1.psf",
            "pdbfile = data/apoa1/apoa1.pdb",
            "rstfile = data/apoa1/apoa1.rst",
        ],
        energy=[
            "switchdist       = 10.0",
            "cutoffdist       = 12.0",
            "pairlistdist     = 13.5",
            "pme_alpha        = 0.34",
            "pme_ngrid_x      = 80",
            "pme_ngrid_y      = 80",
            "pme_ngrid_z      = 64",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=[],
        box=("107.4979", "107.4979", "76.7872"), domain=True, npt_drops_box=True,
    ),
    # ---- CHARMM: UUN ----
    "uun": dict(
        ff="CHARMM",
        input=[
            "topfile = data/uun/toppar/top_all36_prot.rtf,\\",
            "data/uun/toppar/top_all36_na.rtf,\\",
            "data/uun/toppar/top_all36_lipid.rtf,\\",
            "data/uun/toppar/top_all36_cgenff.rtf,\\",
            "data/uun/toppar/top_all36_carb.rtf",
            "parfile = data/uun/toppar/par_all36_prot.prm,\\",
            "data/uun/toppar/par_all36_na.prm,\\",
            "data/uun/toppar/par_all36_carb.prm,\\",
            "data/uun/toppar/par_all36_lipid.prm,\\",
            "data/uun/toppar/par_all36_cgenff.prm",
            "strfile = data/uun/toppar/toppar_water_ions.genesis.str",
            "psffile = data/uun/uun.psf",
            "pdbfile = data/uun/uun.pdb",
            "rstfile = data/uun/uun.rst",
        ],
        energy=[
            "switchdist       = 10",
            "cutoffdist       = 12",
            "pairlistdist     = 13.5",
            "pme_ngrid_x      = 128",
            "pme_ngrid_y      = 128",
            "pme_ngrid_z      = 128",
            "pme_nspline      = 4",
        ],
        energy_water="NONE",
        constraints_extra=[],
        box=("126.5795", "126.5795", "130.6978"), domain=True, npt_drops_box=False,
    ),
    # ---- AMBER: FactorIX ----
    "factorix": dict(
        ff="AMBER",
        input=[
            "prmtopfile = data/factorix/FactorIX.prmtop",
            "ambcrdfile = data/factorix/FactorIX.inpcrd",
            "rstfile    = data/factorix/rst",
        ],
        energy=[
            "switchdist       = 8",
            "cutoffdist       = 8",
            "pairlistdist     = 9.5",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=["water_model      = WAT"],
        box=("142.0855468", "83.3368905", "78.6783548"), domain=False, npt_drops_box=True,
        hmr_topology=True,
    ),
    # ---- GROAMBER: BPTI ----
    "bpti": dict(
        ff="GROAMBER",
        input=[
            "grotopfile = data/bpti/bpti.top",
            "grocrdfile = data/bpti/bpti.gro",
            "rstfile    = data/bpti/rst",
        ],
        energy=[
            "switchdist       = 12.0",
            "cutoffdist       = 12.0",
            "pairlistdist     = 14.0",
            "dielec_const     = 1.0",
            "pme_alpha        = 0.34",
            "pme_ngrid_x      = 64",
            "pme_ngrid_y      = 64",
            "pme_ngrid_z      = 64",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=["water_model      = SOL"],
        box=None, domain=True, npt_drops_box=False,
    ),
    # ---- CHARMM: DPPC ----
    "dppc": dict(
        ff="CHARMM",
        input=[
            "topfile          = data/dppc/top_all36_lipid.rtf",
            "parfile          = data/dppc/par_all36_lipid.prm",
            "strfile          = data/dppc/toppar_water_ions.str",
            "rstfile          = data/dppc/rst",
            "psffile          = data/dppc/dppc.psf",
            "pdbfile          = data/dppc/dppc.pdb",
        ],
        energy=[
            "switchdist       = 10.0",
            "cutoffdist       = 12.0",
            "pairlistdist     = 13.5",
            "table_density    = 20.0",
            "pme_alpha        = 0.34",
            "pme_ngrid_x      = 72",
            "pme_ngrid_y      = 72",
            "pme_ngrid_z      = 72",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=[],
        box=None, domain=True, npt_drops_box=False,
    ),
    # ---- AMBER: AKE ----
    "ake": dict(
        ff="AMBER",
        input=[
            "prmtopfile       = data/ake/ake.top",
            "ambcrdfile       = data/ake/ake.rst",
            "ambreffile       = data/ake/ake.crd",
            "rstfile          = data/ake/restart.rst",
        ],
        energy=[
            "switchdist       = 10.0",
            "cutoffdist       = 10.0",
            "pairlistdist     = 12.0",
            "table_density    = 20.0",
            "pme_alpha        = 0.34",
            "pme_ngrid_x      = 72",
            "pme_ngrid_y      = 72",
            "pme_ngrid_z      = 72",
            "pme_nspline      = 4",
        ],
        energy_water="WAT",
        constraints_extra=["water_model      = WAT"],
        box=None, domain=True, npt_drops_box=False,
    ),
    # ---- AMBER: STMV (huge) ----
    "stmv": dict(
        ff="AMBER",
        input=[
            "prmtopfile = data/stmv/prmtop",
            "ambcrdfile = data/stmv/inpcrd",
        ],
        energy=[
            "switchdist       = 8",
            "cutoffdist       = 8",
            "pairlistdist     = 9.5",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=["water_model      = WAT"],
        box=("221.1723142", "223.1988809", "224.4925841"), domain=False, npt_drops_box=False,
        hmr_topology=True,
    ),
}

ENSEMBLES = ("nve", "nvt", "npt")
DTS = {"2fs": "0.002", "4fs": "0.004"}

# Default steady-state window for a standalone run. The driver overrides nsteps
# and can override eneout_period.
BASE_NSTEPS = 100000
BASE_ENEOUT = 1000
DEFAULT_NBUPDATE_PERIOD = 10
COUPLING_PERIOD_BY_DT = {"2fs": 10, "4fs": 5}


def boundary_block(cfg, ens):
    lines = ["[BOUNDARY]", "type             = PBC"]
    if cfg["box"] is not None and not (ens == "npt" and cfg["npt_drops_box"]):
        bx, by, bz = cfg["box"]
        lines += ["box_size_x       = %s" % bx,
                  "box_size_y       = %s" % by,
                  "box_size_z       = %s" % bz]
    if cfg["domain"]:
        lines += ["domain_x         = 1", "domain_y         = 1", "domain_z         = 1"]
    return lines


def make_input(sysname, ens, dt):
    cfg = SYSTEMS[sysname]
    coupling_period = COUPLING_PERIOD_BY_DT[dt]
    L = []
    # [INPUT]
    L.append("[INPUT]")
    L += cfg["input"]
    L.append("")
    # [OUTPUT] (empty: no trajectory/restart output -> fast, no extra divisibility checks)
    L.append("[OUTPUT]")
    L.append("")
    # [ENERGY]
    L.append("[ENERGY]")
    L.append("forcefield       = %s" % cfg["ff"])
    L.append("electrostatic    = PME")
    L += cfg["energy"]
    L.append("nonbond_kernel   = GPU")
    if cfg["energy_water"] is not None:
        L.append("water_model      = %s" % cfg["energy_water"])
    L.append("")
    # [DYNAMICS]
    L.append("[DYNAMICS]")
    L.append("integrator       = VVER")
    L.append("iseed            = 314159")
    L.append("nsteps           = %d" % BASE_NSTEPS)
    L.append("timestep         = %s" % DTS[dt])
    if dt == "4fs" and not cfg.get("hmr_topology", False):
        L.append("hydrogen_mr      = YES")
        L.append("hmr_target       = all")
        L.append("hmr_ratio        = 3.0")
    L.append("eneout_period    = %d" % BASE_ENEOUT)
    L.append("nbupdate_period  = %d" % DEFAULT_NBUPDATE_PERIOD)
    if ens in ("nvt", "npt"):
        L.append("thermostat_period = %d" % coupling_period)
    if ens == "npt":
        L.append("barostat_period  = %d" % coupling_period)
        L.append("baroscale_period = %d" % coupling_period)
    L.append("")
    # [CONSTRAINTS]
    L.append("[CONSTRAINTS]")
    L.append("rigid_bond       = YES")
    L += cfg["constraints_extra"]
    if dt == "4fs" or cfg.get("hmr_topology", False):
        L.append("hydrogen_mass_upper_bound = 3.3")
    L.append("")
    # [ENSEMBLE]
    L.append("[ENSEMBLE]")
    L.append("ensemble         = %s" % ens.upper())
    L.append("tpcontrol        = %s" % ("NO" if ens == "nve" else "BUSSI"))
    L.append("temperature      = 300")
    if ens == "npt":
        L.append("pressure         = 1.0")
    L.append("group_tp         = YES")
    L.append("")
    # [BOUNDARY]
    L += boundary_block(cfg, ens)
    L.append("")
    return "\n".join(L) + "\n"


def main():
    os.makedirs(OUT, exist_ok=True)
    n = 0
    for sysname in SYSTEMS:
        for ens in ENSEMBLES:
            for dt in DTS:
                fname = "%s_%s_%s.inp" % (sysname, ens, dt)
                with open(os.path.join(OUT, fname), "w") as f:
                    f.write(make_input(sysname, ens, dt))
                n += 1
    print("wrote %d input files to %s" % (n, OUT))


if __name__ == "__main__":
    main()
