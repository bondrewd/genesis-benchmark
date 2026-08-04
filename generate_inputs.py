#!/usr/bin/env python3
"""Generate the GENESIS benchmark input matrix.

The generated inputs cover 9 systems, 3 ensembles, and 2 time steps. Each input
uses paths relative to the benchmark repository root because the benchmark
driver launches GENESIS from that directory.
"""

import os


HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, "inputs")

SYSTEMS = {
    "dhfr": dict(
        forcefield="AMBER",
        input=[
            "prmtopfile = data/dhfr/prmtop",
            "ambcrdfile = data/dhfr/inpcrd",
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
        box=("126.5795", "126.5795", "130.6978"),
        domain=True,
        npt_drops_box=False,
    ),
    "factorix": dict(
        forcefield="AMBER",
        input=[
            "prmtopfile = data/factorix/FactorIX.prmtop",
            "ambcrdfile = data/factorix/FactorIX.inpcrd",
        ],
        energy=[
            "switchdist       = 8",
            "cutoffdist       = 8",
            "pairlistdist     = 9.5",
            "pme_nspline      = 4",
        ],
        box=("142.0855468", "83.3368905", "78.6783548"),
        domain=False,
        npt_drops_box=False,
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
        box=("221.1723142", "223.1988809", "224.4925841"),
        domain=False,
        npt_drops_box=False,
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
        box=("259.2299548", "124.5580494", "123.5021394"),
        domain=False,
        npt_drops_box=False,
    ),
}

ENSEMBLES = ("nve", "nvt", "npt")
TIME_STEPS = {"2fs": "0.002", "4fs": "0.004"}
BASE_NUM_STEPS = 100000
BASE_ENEOUT_PERIOD = 1000
COUPLING_PERIOD_BY_TIME_STEP = {"2fs": 10, "4fs": 5}


def align_assignments(lines):
    """Return input text with every parameter assignment aligned globally."""
    keys = [line.partition("=")[0].strip() for line in lines if "=" in line]
    width = max(len(key) for key in keys)
    aligned = []
    for line in lines:
        if "=" not in line:
            aligned.append(line)
            continue
        key, _, value = line.partition("=")
        aligned.append("%s = %s" % (key.strip().ljust(width), value.strip()))
    return "\n".join(aligned).rstrip() + "\n"


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

    lines.append("[ENERGY]")
    lines.append("forcefield       = %s" % config["forcefield"])
    lines.append("electrostatic    = PME")
    lines += config["energy"]
    lines.append("")

    lines.append("[DYNAMICS]")
    lines.append("integrator       = VVER")
    lines.append("iseed            = 314159")
    lines.append("nsteps           = %d" % BASE_NUM_STEPS)
    lines.append("timestep         = %s" % TIME_STEPS[time_step])
    if time_step == "4fs":
        lines.append("hydrogen_mr      = YES")
        lines.append("hmr_target       = all")
        lines.append("hmr_ratio        = 3.0")
    lines.append("eneout_period    = %d" % BASE_ENEOUT_PERIOD)
    if ensemble in ("nvt", "npt"):
        lines.append("thermostat_period = %d" % coupling_period)
    if ensemble == "npt":
        lines.append("barostat_period  = %d" % coupling_period)
        lines.append("baroscale_period = %d" % coupling_period)
    lines.append("")

    lines.append("[CONSTRAINTS]")
    lines.append("rigid_bond       = YES")
    lines.append("cons_scheme      = MSHAKE")
    lines.append("iter_solute      = 3")
    lines.append("iter_water       = 3")
    if time_step == "4fs":
        lines.append("hydrogen_mass_upper_bound = 3.3")
    lines.append("")

    lines.append("[ENSEMBLE]")
    lines.append("ensemble         = %s" % ensemble.upper())
    lines.append("tpcontrol        = %s" % ("NO" if ensemble == "nve" else "BUSSI"))
    lines.append("temperature      = 300")
    if ensemble == "npt":
        lines.append("pressure         = 1.0")
    lines.append("")

    lines += boundary_block(config, ensemble)
    lines.append("")
    return align_assignments(lines)


def main():
    """Write all generated input files."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    expected_names = {
        "%s_%s_%s.inp" % (system_name, ensemble, time_step)
        for system_name in SYSTEMS
        for ensemble in ENSEMBLES
        for time_step in TIME_STEPS
    }
    for file_name in os.listdir(OUTPUT_DIR):
        if file_name.endswith(".inp") and file_name not in expected_names:
            os.remove(os.path.join(OUTPUT_DIR, file_name))

    files_written = 0
    for system_name in SYSTEMS:
        for ensemble in ENSEMBLES:
            for time_step in TIME_STEPS:
                file_name = "%s_%s_%s.inp" % (system_name, ensemble, time_step)
                output_path = os.path.join(OUTPUT_DIR, file_name)
                with open(output_path, "w") as input_file:
                    input_file.write(make_input(system_name, ensemble, time_step))
                files_written += 1
    print("wrote %d input files to %s" % (files_written, OUTPUT_DIR))


if __name__ == "__main__":
    main()
