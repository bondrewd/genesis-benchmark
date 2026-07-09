#!/usr/bin/env python3
"""Generate the GENESIS benchmark input matrix.

The generated inputs cover 10 systems, 3 ensembles, and each system's supported
time steps. Each generated input uses paths relative to the benchmark repository
root because the benchmark driver launches GENESIS from that directory.
"""

import os


HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, "inputs")

SYSTEMS = {
    "dhfr_27k": dict(
        forcefield="AMBER",
        input=[
            "prmtopfile = data/dhfr_27k/step3_input.parm7",
            "ambcrdfile = data/dhfr_27k/step3_input.rst7",
            "rstfile    = data/dhfr_27k/equil.rst",
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
        box=("1", "1", "1"),
        domain=False,
        npt_drops_box=False,
    ),
    "dhfr_23k": dict(
        forcefield="AMBER",
        input=[
            "prmtopfile = data/dhfr_23k/prmtop",
            "ambcrdfile = data/dhfr_23k/inpcrd",
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
        box=("62.23", "62.23", "62.23"),
        domain=False,
        npt_drops_box=False,
    ),
    "apoa1": dict(
        forcefield="CHARMM",
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
        box=("107.4979", "107.4979", "76.7872"),
        domain=True,
        npt_drops_box=True,
    ),
    "uun": dict(
        forcefield="CHARMM",
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
        box=("126.5795", "126.5795", "130.6978"),
        domain=True,
        npt_drops_box=False,
    ),
    "factorix": dict(
        forcefield="AMBER",
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
        box=("142.0855468", "83.3368905", "78.6783548"),
        domain=False,
        npt_drops_box=True,
        hmr_topology=True,
        time_steps=("4fs",),
    ),
    "bpti": dict(
        forcefield="GROAMBER",
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
        box=None,
        domain=True,
        npt_drops_box=False,
    ),
    "dppc": dict(
        forcefield="CHARMM",
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
        box=None,
        domain=True,
        npt_drops_box=False,
    ),
    "ake": dict(
        forcefield="AMBER",
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
        box=None,
        domain=True,
        npt_drops_box=False,
    ),
    "stmv": dict(
        forcefield="AMBER",
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
        box=("221.1723142", "223.1988809", "224.4925841"),
        domain=False,
        npt_drops_box=False,
        hmr_topology=True,
        time_steps=("4fs",),
    ),
    "cellulose": dict(
        forcefield="AMBER",
        input=[
            "prmtopfile = data/cellulose/prmtop",
            "ambcrdfile = data/cellulose/inpcrd",
        ],
        energy=[
            "switchdist       = 8",
            "cutoffdist       = 8",
            "pairlistdist     = 9.5",
            "pme_nspline      = 4",
        ],
        energy_water=None,
        constraints_extra=["water_model      = WAT"],
        box=("259.2299548", "124.5580494", "123.5021394"),
        domain=False,
        npt_drops_box=False,
    ),
}

ENSEMBLES = ("nve", "nvt", "npt")
TIME_STEPS = {"2fs": "0.002", "4fs": "0.004"}
BASE_NUM_STEPS = 100000
BASE_ENEOUT_PERIOD = 1000
DEFAULT_NBUPDATE_PERIOD = 10
COUPLING_PERIOD_BY_TIME_STEP = {"2fs": 10, "4fs": 5}


def system_time_steps(config):
    """Return supported time-step labels for one system."""
    return config.get("time_steps", tuple(TIME_STEPS))


def boundary_block(config, ensemble):
    """Return the [BOUNDARY] block for one system and ensemble."""
    lines = ["[BOUNDARY]", "type             = PBC"]
    if config["box"] is not None and not (ensemble == "npt" and config["npt_drops_box"]):
        box_x, box_y, box_z = config["box"]
        lines += [
            "box_size_x       = %s" % box_x,
            "box_size_y       = %s" % box_y,
            "box_size_z       = %s" % box_z,
        ]
    if config["domain"]:
        lines += ["domain_x         = 1", "domain_y         = 1", "domain_z         = 1"]
    return lines


def make_input(system_name, ensemble, time_step):
    """Return the GENESIS input text for one benchmark cell."""
    config = SYSTEMS[system_name]
    coupling_period = COUPLING_PERIOD_BY_TIME_STEP[time_step]
    lines = []

    lines.append("[INPUT]")
    lines += config["input"]
    lines.append("")

    lines.append("[OUTPUT]")
    lines.append("")

    lines.append("[ENERGY]")
    lines.append("forcefield       = %s" % config["forcefield"])
    lines.append("electrostatic    = PME")
    lines += config["energy"]
    lines.append("nonbond_kernel   = GPU")
    if config["energy_water"] is not None:
        lines.append("water_model      = %s" % config["energy_water"])
    lines.append("")

    lines.append("[DYNAMICS]")
    lines.append("integrator       = VVER")
    lines.append("iseed            = 314159")
    lines.append("nsteps           = %d" % BASE_NUM_STEPS)
    lines.append("timestep         = %s" % TIME_STEPS[time_step])
    has_hmr_topology = config.get("hmr_topology", False)
    if time_step == "4fs" and not has_hmr_topology:
        lines.append("hydrogen_mr      = YES")
        lines.append("hmr_target       = all")
        lines.append("hmr_ratio        = 3.0")
    lines.append("eneout_period    = %d" % BASE_ENEOUT_PERIOD)
    lines.append("nbupdate_period  = %d" % DEFAULT_NBUPDATE_PERIOD)
    if ensemble in ("nvt", "npt"):
        lines.append("thermostat_period = %d" % coupling_period)
    if ensemble == "npt":
        lines.append("barostat_period  = %d" % coupling_period)
        lines.append("baroscale_period = %d" % coupling_period)
    lines.append("")

    lines.append("[CONSTRAINTS]")
    lines.append("rigid_bond       = YES")
    lines += config["constraints_extra"]
    if time_step == "4fs" or has_hmr_topology:
        lines.append("hydrogen_mass_upper_bound = 3.3")
    lines.append("")

    lines.append("[ENSEMBLE]")
    lines.append("ensemble         = %s" % ensemble.upper())
    lines.append("tpcontrol        = %s" % ("NO" if ensemble == "nve" else "BUSSI"))
    lines.append("temperature      = 300")
    if ensemble == "npt":
        lines.append("pressure         = 1.0")
    lines.append("group_tp         = YES")
    lines.append("")

    lines += boundary_block(config, ensemble)
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main():
    """Write all generated input files."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    files_written = 0
    for system_name, config in SYSTEMS.items():
        for ensemble in ENSEMBLES:
            for time_step in system_time_steps(config):
                file_name = "%s_%s_%s.inp" % (system_name, ensemble, time_step)
                output_path = os.path.join(OUTPUT_DIR, file_name)
                with open(output_path, "w") as input_file:
                    input_file.write(make_input(system_name, ensemble, time_step))
                files_written += 1
    print("wrote %d input files to %s" % (files_written, OUTPUT_DIR))


if __name__ == "__main__":
    main()
