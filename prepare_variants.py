#!/usr/bin/env python3
"""Prepare source-preserving benchmark variants that require topology edits.

The canonical benchmark matrix has exactly one native/archive physical model
per system.  DHFR and ApoA1 require explicit counterions for PME.  Native C36
systems are prepared when restart-derived coordinates, flattened CHARMM cards,
or a static NAMD HMR topology are required.  Neutral native AMBER archives are
prepared when their restart box must be authoritative or static HMR assets are
needed.  No preparation changes the native force-field or water parameters.

Only transformations with complete cross-format validation are enabled.  A
requested unfinished transformation fails closed instead of being relabelled.
Run with the OpenMM environment, which also contains ParmEd::

    conda run -n OpenMM python prepare_variants.py --systems dhfr
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
from dataclasses import dataclass
from decimal import Decimal
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import sys
import tarfile
import tempfile
import warnings


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DEFAULT_OUTPUT_ROOT = DATA / "variants"
DEFAULT_AMBER_TIP3P_IONS = Path(
    "/home/diego/miniforge3/envs/AmberTools26/dat/leap/parm/frcmod.ionsjc_tip3p"
)
SCRIPT_VERSION = 1
PREPARATION_SCRIPT_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
EXPECTED_OPENMM_VERSION = "8.5.2"
EXPECTED_OPENMM_GIT_REVISION = "36a30cbca54e727b216b606f3c011b67201eb8b4"
EXPECTED_OPENMM_RELEASE = False
EXPECTED_PARMED_VERSION = "4.3.1"
AMBERTOOLS26_TIP3P_ION_SOURCE = (
    "ambertools26-dac-26.0.0:dat/leap/parm/frcmod.ionsjc_tip3p"
)
AMBERTOOLS26_TIP3P_ION_SHA256 = (
    "eb7a1a59bf356581f89534287dd59d1bebbb2e5de7176d7c7dd69cb1d4928a55"
)
FORMAL_CHARGE_TOLERANCE_E = 1.0e-3
POST_REPLACEMENT_TOLERANCE_E = 1.0e-6
CHARMM_EXT_COORDINATE_TOLERANCE_ANGSTROM = 5.1e-11
HMR_HYDROGEN_MASS_DA = 3.024
MIN_ION_SEPARATION_NM = 0.5
MIN_ION_SOLUTE_DISTANCE_NM = 0.2
NEUTRALIZATION_SEED_PREFIX = "genesis-benchmark-neutralize-v1"
WATER_NAMES = frozenset(("HOH", "WAT", "SOL", "TIP3", "TIP3P"))
CHARMM_TERM_FIELDS = {
    "bonds": ("atom1", "atom2"),
    "angles": ("atom1", "atom2", "atom3"),
    "dihedrals": ("atom1", "atom2", "atom3", "atom4"),
    "impropers": ("atom1", "atom2", "atom3", "atom4"),
    "donors": ("atom1", "atom2"),
    "acceptors": ("atom1", "atom2"),
    "cmaps": ("atom1", "atom2", "atom3", "atom4", "atom5"),
}


class PreparationError(RuntimeError):
    """Raised when a requested physical transformation is not proven safe."""


@dataclass(frozen=True)
class NativeVariant:
    system: str
    name: str
    family: str
    forcefield: str
    water_model: str
    status: str
    reason: str = ""


VARIANTS = (
    NativeVariant(
        "dhfr", "native_amber_jac_tip3p", "AMBER", "archive-AMBER-JAC", "TIP3P", "ready",
    ),
    NativeVariant(
        "apoa1", "native_charmm27_mtip3p", "CHARMM", "CHARMM27", "mTIP3P", "ready",
    ),
    NativeVariant("uun", "native_charmm36_mtip3p", "CHARMM", "CHARMM36", "mTIP3P", "ready"),
    NativeVariant(
        "factorix", "native_archive_amber_factorix_3site", "AMBER", "custom-factorix",
        "TIP3P", "ready",
    ),
    NativeVariant("bpti", "native_amber03_tip3p", "AMBER", "Amber03", "TIP3P", "ready"),
    NativeVariant(
        "dppc", "native_charmm36lipid_mtip3p", "CHARMM", "CHARMM36lipid", "mTIP3P", "ready",
    ),
    NativeVariant(
        "ake", "native_archive_amber_ake_3site", "AMBER", "archive-AMBER",
        "TIP3P", "ready",
    ),
    NativeVariant(
        "stmv", "native_archive_amber_stmv_3site", "AMBER", "archive-AMBER",
        "TIP3P", "ready",
    ),
    NativeVariant(
        "cellulose", "native_archive_amber_cellulose_3site", "AMBER", "archive-AMBER",
        "TIP3P", "ready",
    ),
)


# These two formally neutral C36 archives need no physical topology mutation.
# Preparation replaces their low-precision PDB coordinate route with the
# authoritative GENESIS restart, flattens the embedded water/ion stream, and
# creates the named-type and solute-only HMR PSFs required by NAMD.  The pinned
# hashes are the outputs of the independently validated surgical prototype.
NATIVE_CHARMM_ARCHIVES = {
    "dppc": {
        "forcefield": "CHARMM36lipid",
        "psf": "dppc.psf",
        "pdb": "dppc.pdb",
        "restart": "rst",
        "topologies": ("top_all36_lipid.rtf",),
        "parameters": ("par_all36_lipid.prm",),
        "stream": "toppar_water_ions.str",
        "atom_count": 36126,
        "box_angstrom": (69.47812837130373, 69.47812837130373, 71.64969487712301),
        "water_count": 5022,
        "existing_ions": {},
        "bonded_term_counts": {
            "bonds": 35964,
            "angles": 45522,
            "proper_dihedrals": 56538,
            "impropers": 324,
            "cmaps": 0,
        },
        "exception_count": 133002,
        "hmr_solute_hydrogens": 12960,
        "hmr_donors": 6156,
        "normal_total_mass_dalton": "209390.44320",
        "asset_sha256": {
            "par_all36_lipid.prm": "ffd11de3b382738601dd4139d43fe3e51f39c9c31bd2fea54dab89b506db627a",
            "system.coor": "e23b75c05594482646a46797f7449895af3829049822cddc1c5b25fd8bfb450c",
            "system.crd": "ad6ca4407bfc48fb266d2ecdb32eea1581aeabeeaa25f5e0362fb1829d2367d2",
            "system.pdb": "c41b09ebdba9460bd83487bb08a48b2a814aaed1a4c87efc4da862a242850a9f",
            "system.psf": "4e2aae4e420c7d636696b0579f1cb05f1fffb615676c92594fd0f79c8c6db3e3",
            "system_hmr_xplor.psf": "0e89c2fed7f4cecf4638cd504e577235c82ffcad4e943983dcc814327bad9322",
            "system_xplor.psf": "e84623c2540fe809d513b8e7db1f13aca87e0982d653ab7ff51b7ff09cab9221",
            "top_all36_lipid.rtf": "7615982187e2e81f878ec982ad6b56a872543722988c17a12a63df6895224dfb",
            "water_ions.prm": "8d381380689fb0b01c8f1f56f6730a589c2d6d1bd1620b2e61482bbc819689bf",
            "water_ions.rtf": "7e4e9cc33aee29488bcac3aa61e09ec3429415bb33d527bbcf64cd6794fb3f86",
            "water_ions_nbfix.prm": "a3d5825d3f1a7e438f1bcc2465ce0307ce59b5e0ef6ac44d7a787244d92f31f7",
        },
    },
    "uun": {
        "forcefield": "CHARMM36",
        "psf": "uun.psf",
        "pdb": "uun.pdb",
        "restart": "uun.rst",
        "topologies": tuple("toppar/top_all36_%s.rtf" % name for name in (
            "prot", "na", "lipid", "cgenff", "carb",
        )),
        "parameters": tuple("toppar/par_all36_%s.prm" % name for name in (
            "prot", "na", "carb", "lipid", "cgenff",
        )),
        "stream": "toppar/toppar_water_ions.genesis.str",
        "atom_count": 216726,
        "box_angstrom": (126.57952880859375, 126.57952880859375, 130.69783020019531),
        "water_count": 45106,
        "existing_ions": {"CLA": 104, "POT": 192},
        "bonded_term_counts": {
            "bonds": 216258,
            "angles": 198186,
            "proper_dihedrals": 215952,
            "impropers": 4784,
            "cmaps": 1456,
        },
        "exception_count": 584178,
        "hmr_solute_hydrogens": 46984,
        "hmr_donors": 24680,
        "normal_total_mass_dalton": "1317099.00280",
        "asset_sha256": {
            "par_all36_carb.prm": "ff3bc65cf1bb615ee3280140cd813b7d5f1d4838a79b8c89814ce71cc352a7e5",
            "par_all36_cgenff.prm": "420cfb4355cd8bdb3a1e0d03d2b14c903251b201e1fd366c4e46b7d515ae7258",
            "par_all36_lipid.prm": "b0ffe75adbb42695da6cfd2b4b083a91156f048f931cecaeaecf47974eea60ec",
            "par_all36_na.prm": "5ad1e4acd8e95590c6da4400e81f2a9343362394cb040202ddee984a15377f6d",
            "par_all36_prot.prm": "0e6eedd83c604d94de88071461373000c1beb91a8087fec6b6eca38341a86c4a",
            "system.coor": "fd4569f5268db9ba37edfef06861d0b19f3562848404d453812de5be734fba7b",
            "system.crd": "4ed8a40c99ca465e68af6f9bfeaacd8ab5e9f2eb4713e451c58b0dee52adb550",
            "system.pdb": "2ce360fd32b49351af0f28725079cc5830408cd77e61a131fa568ee92b646470",
            "system.psf": "5fa25fc8e9eaa87abad88d08e467bb440df924dc9d38c90e32b98fb8c99afbef",
            "system_hmr_xplor.psf": "80e3b64682a1b622dafbd5438b06669a6458819bffc35cc8f0b125c7fd1da919",
            "system_xplor.psf": "5fa25fc8e9eaa87abad88d08e467bb440df924dc9d38c90e32b98fb8c99afbef",
            "top_all36_carb.rtf": "f92a53e3f1e9718822e3a0c45d2d7df953180ce4de7fb21b5854a83f6aca43c1",
            "top_all36_cgenff.rtf": "3fd582686faea5aefd0f1f64f7bc9ea07adc2afa592ebedcb30bf15d099cad01",
            "top_all36_lipid.rtf": "df667f92758cd5a3eb0596e925c4b3133ac09b36904619fd376c4b5ff02b797d",
            "top_all36_na.rtf": "42c8e5621f5f1fd5b75df14de02fb07fbc16fdbed76ebce4bf40cbc436310e24",
            "top_all36_prot.rtf": "cf151676a676c744d1fc0267851244a9a35abebc297e582357550b11eeceedd1",
            "water_ions.prm": "228c8307c7528d22f4a26fd8ebfa0aa1ecce0a374d44623ea3565aa724ba5b74",
            "water_ions.rtf": "7e4e9cc33aee29488bcac3aa61e09ec3429415bb33d527bbcf64cd6794fb3f86",
            "water_ions_nbfix.prm": "c25b5274862af51d44cdb7bc45bc974939b54f68d4bb69368b3117697c2a9798",
        },
    },
}


# Formally neutral native AMBER systems need no force-field mutation.  Their
# preparation copies the authoritative restart verbatim, reconciles only the
# prmtop BOX_DIMENSIONS field when it is stale, and creates solute-only static
# HMR topologies for engines that cannot repartition masses at runtime.
# GROMACS is deliberately opt-in per system: an asset is advertised only after
# an actual-written-file Amber/GROMACS single point meets the strict energy and
# force thresholds in numerical_amber_gromacs_validation.
NATIVE_AMBER_ARCHIVES = {
    "bpti": {
        "topology": "prmtop",
        "coordinates": "inpcrd",
        "atom_count": 27712,
        "box_angstrom": (65.3318, 65.3318, 65.3318),
        "gromacs_equivalence": True,
    },
    "ake": {
        "topology": "ake.top",
        "coordinates": "ake.rst",
        "atom_count": 62475,
        "box_angstrom": (86.4941345, 83.2631699, 87.1693322),
        "gromacs_equivalence": True,
    },
    "factorix": {
        "topology": "FactorIX.prmtop",
        "coordinates": "FactorIX.inpcrd",
        "atom_count": 90906,
        "box_angstrom": (142.0855468, 83.3368905, 78.6783548),
        "gromacs_equivalence": True,
    },
    "cellulose": {
        "topology": "prmtop",
        "coordinates": "inpcrd",
        "atom_count": 408609,
        "box_angstrom": (259.2299548, 124.5580494, 123.5021394),
        "gromacs_equivalence": True,
    },
    "stmv": {
        "topology": "prmtop",
        "coordinates": "inpcrd",
        "atom_count": 1067095,
        "box_angstrom": (221.1723142, 223.1988809, 224.4925841),
        "gromacs_equivalence": True,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json(payload: object) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n"


def integer_charge(charge: float, label: str) -> int:
    formal = int(math.floor(charge + 0.5))
    if not math.isfinite(charge) or abs(charge - formal) > FORMAL_CHARGE_TOLERANCE_E:
        raise PreparationError(
            "%s has nonintegral charge %.12g e (tolerance %.1g e)"
            % (label, charge, FORMAL_CHARGE_TOLERANCE_E)
        )
    return formal


def parse_csv(values: list[str] | None) -> set[str] | None:
    if not values:
        return None
    result = set()
    for value in values:
        result.update(item.strip() for item in value.split(",") if item.strip())
    return result


def unpack_archive(system: str, parent: Path) -> Path:
    archive = DATA / (system + ".tgz")
    if not archive.is_file():
        raise PreparationError("missing source archive: %s" % archive)
    with tarfile.open(archive, "r:gz") as bundle:
        members = bundle.getmembers()
        root = parent.resolve()
        for member in members:
            target = (parent / member.name).resolve()
            if target != root and root not in target.parents:
                raise PreparationError("unsafe member %r in %s" % (member.name, archive))
            if member.issym() or member.islnk():
                raise PreparationError("links are not allowed in %s" % archive)
        bundle.extractall(parent, members=members, filter="data")
    source = parent / system
    if not source.is_dir():
        raise PreparationError("%s does not contain %s/" % (archive, system))
    return source


def is_water(residue: object) -> bool:
    atoms = list(residue.atoms())
    symbols = sorted(atom.element.symbol if atom.element is not None else "EP" for atom in atoms)
    return residue.name.upper() in WATER_NAMES and symbols == ["H", "H", "O"]


def is_monatomic_ion(residue: object) -> bool:
    atoms = list(residue.atoms())
    return (
        len(atoms) == 1 and atoms[0].element is not None
        and atoms[0].element.symbol in ("Na", "Cl", "K", "Ca", "Mg")
    )


def orthorhombic_lengths_nm(topology: object) -> tuple[float, float, float]:
    from openmm import unit

    vectors = topology.getPeriodicBoxVectors()
    if vectors is None:
        raise PreparationError("explicit neutralization requires periodic box vectors")
    rows = vectors.value_in_unit(unit.nanometer)
    off_diagonal = (rows[0][1], rows[0][2], rows[1][0], rows[1][2], rows[2][0], rows[2][1])
    if any(abs(value) > 1.0e-8 for value in off_diagonal):
        raise PreparationError("neutralization currently requires an orthorhombic box")
    lengths = (float(rows[0][0]), float(rows[1][1]), float(rows[2][2]))
    if any(value <= 0.0 for value in lengths):
        raise PreparationError("periodic box lengths must be positive")
    return lengths


def pbc_distance_squared(a: object, b: object, lengths: tuple[float, float, float]) -> float:
    total = 0.0
    for dimension, length in enumerate(lengths):
        delta = float(a[dimension]) - float(b[dimension])
        delta -= length * math.floor(delta / length + 0.5)
        total += delta * delta
    return total


def select_water_replacements(
    topology: object,
    positions: object,
    system_name: str,
    count: int,
) -> tuple[list[tuple], dict]:
    """Select deterministic, spatially valid water oxygen sites for ions."""
    from openmm import unit
    from openmm.app import element

    coordinates = positions.value_in_unit(unit.nanometer)
    lengths = orthorhombic_lengths_nm(topology)
    seed = "%s:%s" % (NEUTRALIZATION_SEED_PREFIX, system_name)
    waters = [residue for residue in topology.residues() if is_water(residue)]
    candidates = []
    for residue in waters:
        oxygen = next(atom for atom in residue.atoms() if atom.element == element.oxygen)
        # The reproducibility contract hashes the zero-based residue index in
        # the untouched source topology, not a water ordinal or converted index.
        rank = hashlib.sha256(
            seed.encode("utf-8") + b"\0" + str(residue.index).encode("ascii")
        ).digest()
        candidates.append((rank, residue.index, residue, oxygen, coordinates[oxygen.index]))
    candidates.sort(key=lambda item: (item[0], item[1]))

    existing_ions = []
    solute_heavy = []
    existing_ion_counts = {}
    for residue in topology.residues():
        atoms = list(residue.atoms())
        if is_water(residue):
            continue
        if is_monatomic_ion(residue):
            existing_ions.append(coordinates[atoms[0].index])
            existing_ion_counts[residue.name] = existing_ion_counts.get(residue.name, 0) + 1
            continue
        solute_heavy.extend(
            coordinates[atom.index] for atom in atoms
            if atom.element is not None and atom.element != element.hydrogen
        )

    selected = []
    ion_cutoff2 = MIN_ION_SEPARATION_NM * MIN_ION_SEPARATION_NM
    solute_cutoff2 = MIN_ION_SOLUTE_DISTANCE_NM * MIN_ION_SOLUTE_DISTANCE_NM
    for candidate in candidates:
        position = candidate[4]
        if any(pbc_distance_squared(position, other, lengths) < ion_cutoff2
               for other in existing_ions):
            continue
        if any(pbc_distance_squared(position, other[4], lengths) < ion_cutoff2
               for other in selected):
            continue
        if any(pbc_distance_squared(position, other, lengths) < solute_cutoff2
               for other in solute_heavy):
            continue
        selected.append(candidate)
        if len(selected) == count:
            break
    if len(selected) != count:
        raise PreparationError("found only %d of %d safe ion sites" % (len(selected), count))

    metadata = {
        "method": "replace_waters",
        "seed": seed,
        "tolerance_e": FORMAL_CHARGE_TOLERANCE_E,
        "min_ion_separation_nm": MIN_ION_SEPARATION_NM,
        "min_ion_solute_distance_nm": MIN_ION_SOLUTE_DISTANCE_NM,
        "neutralization_ion": "Na+",
        "neutralization_ion_charge_e": 1,
        "neutralization_ion_count": count,
        "replaced_water_count": count,
        "existing_ions": dict(sorted(existing_ion_counts.items())),
        "water_count_before": len(waters),
        "water_count_after": len(waters) - count,
        "selected_water_residue_indices": [item[1] for item in selected],
        "selected_water_oxygen_positions_nm": [list(map(float, item[4])) for item in selected],
    }
    return selected, metadata


def parse_amber_na_parameters(path: Path) -> tuple[float, float, float]:
    """Read Na+ mass, Rmin/2, and epsilon from an Amber frcmod file."""
    if not path.is_file():
        raise PreparationError("missing Amber TIP3P ion parameters: %s" % path)
    section = None
    mass = radius = epsilon = None
    for raw_line in path.read_text(encoding="ascii").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if line in ("MASS", "NONBON"):
            section = line
            continue
        if not line or not line.startswith("Na+"):
            continue
        fields = line.split()
        if section == "MASS" and len(fields) >= 2:
            mass = float(fields[1])
        elif section == "NONBON" and len(fields) >= 3:
            radius, epsilon = float(fields[1]), float(fields[2])
    if mass is None or radius is None or epsilon is None:
        raise PreparationError("could not read Na+ MASS/NONBON entries from %s" % path)
    return mass, radius, epsilon


def term_connectivity(
    structure: object, selected_atom_indices: set[int],
) -> dict[str, Counter[tuple[int, ...]]]:
    def retained(items: object, names: tuple[str, ...]) -> Counter[tuple[int, ...]]:
        result = Counter()
        for item in items:
            indices = tuple(getattr(item, name).idx for name in names)
            if not any(index in selected_atom_indices for index in indices):
                result[indices] += 1
        return result

    return {
        "bonds": retained(structure.bonds, ("atom1", "atom2")),
        "angles": retained(structure.angles, ("atom1", "atom2", "atom3")),
        "dihedrals": retained(structure.dihedrals, ("atom1", "atom2", "atom3", "atom4")),
    }


def mutate_native_dhfr(
    source_topology: Path,
    source_coordinates: Path,
    ion_parameters: Path,
) -> tuple[object, object, dict, dict]:
    """Replace selected native TIP3P waters in an AmberParm with JC Na+."""
    import parmed as pmd
    from parmed.amber.mask import AmberMask
    from parmed.topologyobjects import AtomType
    from parmed.tools.actions import addLJType
    from parmed.tools.argumentlist import ArgumentList
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile

    source = pmd.load_file(str(source_topology), xyz=str(source_coordinates))
    source_charge = math.fsum(float(atom.charge) for atom in source.atoms)
    source_formal = integer_charge(source_charge, "native DHFR source")
    if source_formal != -11:
        raise PreparationError("native DHFR formal charge changed: expected -11, got %d" % source_formal)

    app_topology = AmberPrmtopFile(str(source_topology)).topology
    app_coordinates = AmberInpcrdFile(str(source_coordinates))
    if app_coordinates.boxVectors is None:
        raise PreparationError("native DHFR coordinate file has no periodic box")
    app_topology.setPeriodicBoxVectors(app_coordinates.boxVectors)
    coordinate_box_angstrom = [
        10.0 * length for length in orthorhombic_lengths_nm(app_topology)
    ]
    source.box = coordinate_box_angstrom + [90.0, 90.0, 90.0]
    if coordinate_delta(
        source.coordinates, app_coordinates.positions.value_in_unit(unit.angstrom)
    ) > 1.0e-6:
        raise PreparationError("ParmEd and OpenMM read different native DHFR coordinates")
    selected, metadata = select_water_replacements(
        app_topology, app_coordinates.positions, "dhfr", -source_formal
    )
    selected_residues = {item[1] for item in selected}
    if max(selected_residues) >= len(source.residues):
        raise PreparationError("OpenMM/ParmEd residue indexing differs for native DHFR")

    native = copy.deepcopy(source)
    sodium_mass, sodium_radius, sodium_epsilon = parse_amber_na_parameters(ion_parameters)
    oxygen_indices = []
    hydrogen_indices = []
    selected_source_atoms = set()
    for residue_index in sorted(selected_residues):
        residue = native.residues[residue_index]
        if residue.name.upper() not in WATER_NAMES or len(residue.atoms) != 3:
            raise PreparationError("selected DHFR residue %d is not native TIP3P water" % residue_index)
        oxygen = [atom for atom in residue.atoms if atom.atomic_number == 8]
        hydrogens = [atom for atom in residue.atoms if atom.atomic_number == 1]
        if len(oxygen) != 1 or len(hydrogens) != 2:
            raise PreparationError("selected DHFR water %d has invalid elements" % residue_index)
        oxygen_indices.append(oxygen[0].idx + 1)
        hydrogen_indices.extend(atom.idx + 1 for atom in hydrogens)
        selected_source_atoms.update(atom.idx for atom in residue.atoms)

    original_lj_radius = tuple(float(value) for value in native.LJ_radius)
    original_lj_depth = tuple(float(value) for value in native.LJ_depth)
    original_terms = term_connectivity(native, selected_source_atoms)
    add_mask = "@" + ",".join(map(str, oxygen_indices))
    addLJType(
        native,
        ArgumentList(
            "%s radius %.12g epsilon %.12g" % (add_mask, sodium_radius, sodium_epsilon)
        ),
    ).execute()

    for residue_index in sorted(selected_residues):
        residue = native.residues[residue_index]
        oxygen = next(atom for atom in residue.atoms if atom.atomic_number == 8)
        residue.name = "Na+"
        oxygen.name = "Na+"
        oxygen.type = "Na+"
        oxygen.charge = 1.0
        oxygen.mass = sodium_mass
        oxygen.atomic_number = 11
    native.strip(AmberMask(native, "@" + ",".join(map(str, hydrogen_indices))))
    native.remake_parm()
    sodium_atom_type = AtomType("Na+", None, sodium_mass, atomic_number=11, charge=1.0)
    sodium_atom_type.set_lj_params(
        sodium_epsilon, sodium_radius, sodium_epsilon, sodium_radius
    )
    for atom in native.atoms:
        if atom.residue.name == "Na+":
            atom.atom_type = sodium_atom_type

    selected_positions_by_residue = {item[1]: item[4] for item in selected}
    for residue_index in selected_residues:
        residue = native.residues[residue_index]
        if len(residue.atoms) != 1 or residue.atoms[0].atomic_number != 11:
            raise PreparationError("native water replacement changed residue ordering")
        observed_nm = native.coordinates[residue.atoms[0].idx] / 10.0
        expected_nm = selected_positions_by_residue[residue_index]
        if any(abs(float(observed_nm[dimension]) - float(expected_nm[dimension])) > 1.0e-8
               for dimension in range(3)):
            raise PreparationError("neutralizing ion moved from its selected water oxygen site")
    remaining_waters = sum(
        residue.name.upper() in WATER_NAMES and len(residue.atoms) == 3
        for residue in native.residues
    )
    if remaining_waters != metadata["water_count_after"]:
        raise PreparationError("native water count changed unexpectedly during neutralization")

    if tuple(float(value) for value in native.LJ_radius[:len(original_lj_radius)]) != original_lj_radius:
        raise PreparationError("adding Na+ changed a pre-existing native LJ radius")
    if tuple(float(value) for value in native.LJ_depth[:len(original_lj_depth)]) != original_lj_depth:
        raise PreparationError("adding Na+ changed a pre-existing native LJ well depth")
    post_charge = math.fsum(float(atom.charge) for atom in native.atoms)
    if integer_charge(post_charge, "neutral native DHFR") != 0 or abs(post_charge) > POST_REPLACEMENT_TOLERANCE_E:
        raise PreparationError("neutralized native DHFR charge is %.12g e" % post_charge)

    metadata.update({
        "source_charge_e": source_charge,
        "source_formal_charge_e": source_formal,
        "pre_neutralization_charge_e": source_charge,
        "pre_neutralization_formal_charge_e": source_formal,
        "post_neutralization_charge_e": post_charge,
        "post_neutralization_formal_charge_e": 0,
    })
    preservation = {
        "original_atom_count": len(source.atoms),
        "neutral_atom_count": len(native.atoms),
        "removed_hydrogen_count": len(hydrogen_indices),
        "original_terms": original_terms,
        "source": source,
        "selected_source_atoms": selected_source_atoms,
        "sodium_mass": sodium_mass,
        "sodium_radius": sodium_radius,
        "sodium_epsilon": sodium_epsilon,
    }
    return native, app_topology, metadata, preservation


def copy_with_hmr(structure: object) -> tuple[object, dict]:
    result = copy.deepcopy(structure)
    changed_hydrogens = set()
    changed_heavy_atoms = set()
    before = math.fsum(float(atom.mass) for atom in result.atoms)
    for atom in result.atoms:
        if atom.atomic_number != 1 or atom.residue.name.upper() in WATER_NAMES:
            continue
        heavy = []
        for bond in atom.bonds:
            partner = bond.atom2 if bond.atom1 is atom else bond.atom1
            if partner.atomic_number != 1:
                heavy.append(partner)
        if len(heavy) != 1:
            raise PreparationError("hydrogen %s has %d bonded heavy atoms" % (atom, len(heavy)))
        delta = HMR_HYDROGEN_MASS_DA - float(atom.mass)
        if delta < -1.0e-6 or float(heavy[0].mass) - delta <= 0.0:
            raise PreparationError("invalid native HMR mass transfer at %s" % atom)
        atom.mass = HMR_HYDROGEN_MASS_DA
        heavy[0].mass -= delta
        changed_hydrogens.add(atom.idx)
        changed_heavy_atoms.add(heavy[0].idx)
    # The legacy source predates the optional ATOMIC_NUMBER flag.  It becomes
    # essential after HMR because element inference from the modified masses is
    # no longer valid (for example a repartitioned nitrogen can resemble Li).
    if "ATOMIC_NUMBER" not in result.parm_data:
        result.add_flag(
            "ATOMIC_NUMBER", "10I8",
            data=[atom.atomic_number for atom in result.atoms], after="CHARGE",
        )
    result.remake_parm()
    after = math.fsum(float(atom.mass) for atom in result.atoms)
    if not changed_hydrogens or abs(after - before) > 1.0e-6:
        raise PreparationError("native HMR failed mass conservation")
    return result, {
        "method": "native-topology heavy-atom mass transfer; rigid waters excluded",
        "hydrogen_mass_dalton": HMR_HYDROGEN_MASS_DA,
        "particles_with_changed_mass": len(changed_hydrogens | changed_heavy_atoms),
        "water_particles_with_changed_mass": 0,
        "total_mass_delta_dalton": after - before,
    }


def normalize_amber_header(path: Path) -> None:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or not lines[0].startswith("%VERSION"):
        raise PreparationError("unexpected Amber topology header in %s" % path)
    lines[0] = "%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/70  00:00:00"
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def normalize_gromacs_header(path: Path) -> None:
    text = path.read_text(encoding="ascii")
    marker = "[ defaults ]"
    if marker not in text:
        raise PreparationError("unexpected GROMACS topology header in %s" % path)
    path.write_text("; Generated by prepare_variants.py\n\n" + marker + text.split(marker, 1)[1], encoding="ascii")


def write_g96(path: Path, structure: object) -> None:
    """Write a high-precision, fixed-format GROMACS-96 coordinate file."""
    coordinates = structure.coordinates / 10.0
    if structure.box is None or len(structure.box) < 6:
        raise PreparationError("G96 export requires a periodic box")
    if any(abs(float(angle) - 90.0) > 1.0e-8 for angle in structure.box[3:6]):
        raise PreparationError("G96 export currently requires an orthorhombic box")
    lines = ["TITLE", "Generated by prepare_variants.py", "END", "POSITION"]
    for atom_number, atom in enumerate(structure.atoms, 1):
        x, y, z = coordinates[atom.idx]
        lines.append(
            "%5d %-5.5s %-5.5s%7d%15.9f%15.9f%15.9f"
            % ((atom.residue.idx + 1) % 100000, atom.residue.name, atom.name,
               atom_number % 10000000, x, y, z)
        )
    lines.extend(("END", "BOX", "%15.9f%15.9f%15.9f" % tuple(
        float(value) / 10.0 for value in structure.box[:3]
    ), "END"))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def read_g96(path: Path) -> tuple[list[tuple[str, str]], list[list[float]], list[float]]:
    """Strictly reload the single-frame G96 subset emitted by :func:`write_g96`.

    Coordinates and box lengths are returned in Angstrom.  Parsing the actual
    serialized asset prevents an in-memory coordinate array from masking file
    precision or ordering errors during cross-engine validation.
    """
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) < 9 or lines[:4] != [
        "TITLE", "Generated by prepare_variants.py", "END", "POSITION",
    ]:
        raise PreparationError("unexpected G96 header in %s" % path)
    try:
        position_end = lines.index("END", 4)
    except ValueError as error:
        raise PreparationError("unterminated G96 POSITION block in %s" % path) from error
    if position_end + 3 >= len(lines) or lines[position_end + 1] != "BOX" or lines[-1] != "END":
        raise PreparationError("unexpected G96 block order in %s" % path)
    if position_end + 4 != len(lines):
        raise PreparationError("unexpected trailing G96 data in %s" % path)

    names = []
    coordinates = []
    for expected_index, line in enumerate(lines[4:position_end], 1):
        if len(line) != 69 or line[5] != " " or line[11] != " ":
            raise PreparationError("invalid fixed-width G96 atom line %d" % expected_index)
        try:
            residue_number = int(line[0:5])
            atom_number = int(line[17:24])
            xyz_nm = [float(line[start:start + 15]) for start in (24, 39, 54)]
        except ValueError as error:
            raise PreparationError("invalid G96 atom line %d" % expected_index) from error
        if residue_number < 0 or atom_number != expected_index % 10000000:
            raise PreparationError("invalid G96 numbering at atom %d" % expected_index)
        names.append((line[6:11].strip(), line[12:17].strip()))
        coordinates.append([10.0 * value for value in xyz_nm])

    box_line = lines[position_end + 2]
    if len(box_line) != 45:
        raise PreparationError("invalid fixed-width G96 box in %s" % path)
    try:
        box = [10.0 * float(box_line[start:start + 15]) for start in (0, 15, 30)]
    except ValueError as error:
        raise PreparationError("invalid G96 box in %s" % path) from error
    if any(value <= 0.0 for value in box):
        raise PreparationError("nonpositive G96 box length in %s" % path)
    return names, coordinates, box


def system_charge(system: object) -> float:
    from openmm import NonbondedForce, unit

    forces = [force for force in system.getForces() if isinstance(force, NonbondedForce)]
    if len(forces) != 1:
        raise PreparationError("expected one NonbondedForce, found %d" % len(forces))
    return math.fsum(
        forces[0].getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
        for index in range(forces[0].getNumParticles())
    )


def mass_structure(structure: object) -> float:
    return math.fsum(float(atom.mass) for atom in structure.atoms)


def mass_system(system: object) -> float:
    from openmm import unit

    return math.fsum(
        system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(system.getNumParticles())
    )


def system_box_angstrom(system: object) -> list[float]:
    from openmm import unit

    vectors = system.getDefaultPeriodicBoxVectors()
    rows = [vector.value_in_unit(unit.angstrom) for vector in vectors]
    if any(abs(float(rows[row][column])) > 1.0e-8
           for row in range(3) for column in range(3) if row != column):
        raise PreparationError("reloaded System has a non-orthorhombic box")
    result = [float(rows[index][index]) for index in range(3)]
    if any(value <= 0.0 for value in result):
        raise PreparationError("reloaded System has a nonpositive box length")
    return result


def system_topology_signature(system: object) -> dict:
    """Return exact index-level topology features that converters must preserve."""
    from openmm import (
        CMAPTorsionForce, HarmonicAngleForce, HarmonicBondForce,
        NonbondedForce, PeriodicTorsionForce, unit,
    )

    counts = {"bonds": 0, "angles": 0, "torsions": 0, "cmaps": 0}
    zero_energy_torsions = 0
    exception_pairs = None
    for force in system.getForces():
        if isinstance(force, HarmonicBondForce):
            counts["bonds"] += force.getNumBonds()
        elif isinstance(force, HarmonicAngleForce):
            counts["angles"] += force.getNumAngles()
        elif isinstance(force, PeriodicTorsionForce):
            for index in range(force.getNumTorsions()):
                force_constant = force.getTorsionParameters(index)[6]
                if abs(force_constant.value_in_unit(unit.kilojoule_per_mole)) <= 1.0e-15:
                    zero_energy_torsions += 1
                else:
                    counts["torsions"] += 1
        elif isinstance(force, CMAPTorsionForce):
            counts["cmaps"] += force.getNumTorsions()
        elif isinstance(force, NonbondedForce):
            if exception_pairs is not None:
                raise PreparationError("reloaded System has multiple NonbondedForce objects")
            exception_pairs = {
                tuple(sorted(force.getExceptionParameters(index)[:2]))
                for index in range(force.getNumExceptions())
            }
    if exception_pairs is None:
        raise PreparationError("reloaded System has no NonbondedForce")
    constraints = {}
    for index in range(system.getNumConstraints()):
        atom1, atom2, distance = system.getConstraintParameters(index)
        pair = tuple(sorted((atom1, atom2)))
        constraints[pair] = distance.value_in_unit(unit.nanometer)
    virtual_sites = tuple(
        index for index in range(system.getNumParticles()) if system.isVirtualSite(index)
    )
    return {
        "bonded_force_term_counts": counts,
        "zero_energy_torsions": zero_energy_torsions,
        "exception_pairs": exception_pairs,
        "constraints": constraints,
        "virtual_sites": virtual_sites,
    }


def coordinate_delta(reference: object, candidate: object) -> float:
    if len(reference) != len(candidate):
        raise PreparationError("coordinate atom counts differ")
    return max(
        abs(float(expected[dimension]) - float(observed[dimension]))
        for expected, observed in zip(reference, candidate) for dimension in range(3)
    )


def asset_record(
    atom_count: int,
    raw_charge: float,
    box: object,
    coordinate_max_delta: float,
    normal_mass: float,
    hmr_mass: float | None,
) -> dict:
    return {
        "status": "pass",
        "atom_count": atom_count,
        "raw_charge_e": raw_charge,
        "formal_charge_e": integer_charge(raw_charge, "reloaded asset"),
        "box_angstrom": [float(value) for value in box[:3]],
        "coordinate_max_delta_angstrom": coordinate_max_delta,
        "normal_total_mass_dalton": normal_mass,
        "hmr_total_mass_dalton": hmr_mass,
    }


def verify_source_preservation(neutral: object, preservation: dict) -> None:
    source = preservation["source"]
    selected_atoms = preservation["selected_source_atoms"]
    expected = []
    for atom in source.atoms:
        if atom.idx not in selected_atoms:
            expected.append(atom)
        elif atom.atomic_number == 8:
            expected.append(None)
    if len(expected) != len(neutral.atoms):
        raise PreparationError("native DHFR atom deletion did not match selected waters")
    for output, original in zip(neutral.atoms, expected):
        if original is None:
            if output.name != "Na+" or output.residue.name != "Na+" or output.atomic_number != 11:
                raise PreparationError("selected native water was not converted to Na+")
            continue
        fields = (
            output.name, output.residue.name, output.type, output.atomic_number, output.nb_idx,
        )
        reference = (
            original.name, original.residue.name, original.type,
            original.atomic_number, original.nb_idx,
        )
        if (fields != reference
                or abs(float(output.charge) - float(original.charge)) > 1.0e-12
                or abs(float(output.mass) - float(original.mass)) > 1.0e-10):
            raise PreparationError("native atom parameters changed outside selected waters")

    output_to_source = {}
    for output, original in zip(neutral.atoms, expected):
        if original is not None:
            output_to_source[output.idx] = original.idx
    selected_output = set(range(len(neutral.atoms))) - set(output_to_source)
    output_terms = term_connectivity(neutral, selected_output)
    remapped = {}
    for name, terms in output_terms.items():
        remapped[name] = Counter()
        for term, count in terms.items():
            remapped[name][tuple(output_to_source[index] for index in term)] += count
    if remapped != preservation["original_terms"]:
        raise PreparationError("bonded connectivity changed outside selected native waters")


def export_native_dhfr(
    directory: Path,
    neutral: object,
    hmr: object,
) -> tuple[dict, object, object]:
    import parmed as pmd
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME

    neutral.save(str(directory / "system.prmtop"), overwrite=True)
    neutral.save(str(directory / "system.inpcrd"), overwrite=True)
    hmr.save(str(directory / "system_hmr.prmtop"), overwrite=True)
    normalize_amber_header(directory / "system.prmtop")
    normalize_amber_header(directory / "system_hmr.prmtop")

    canonical_top = AmberPrmtopFile(str(directory / "system.prmtop"))
    canonical_coordinates = AmberInpcrdFile(str(directory / "system.inpcrd"))
    if canonical_coordinates.boxVectors is None:
        raise PreparationError("exported native DHFR lost its periodic box")
    canonical_top.topology.setPeriodicBoxVectors(canonical_coordinates.boxVectors)
    normal_system = canonical_top.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds, rigidWater=True, removeCMMotion=False,
        ewaldErrorTolerance=1.0e-5,
    )
    runtime_hmr_system = canonical_top.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds, rigidWater=True, removeCMMotion=False,
        ewaldErrorTolerance=1.0e-5,
        hydrogenMass=HMR_HYDROGEN_MASS_DA * unit.dalton,
    )

    write_native_amber_gromacs_topology(neutral, directory / "system.top")
    write_native_amber_gromacs_topology(hmr, directory / "system_hmr.top")
    write_g96(directory / "system.g96", neutral)

    assets = {
        "GENESIS": {
            "format": "AMBER", "topology": "system.prmtop", "coordinates": "system.inpcrd",
            "topology_definitions": [], "parameters": [],
        },
        "OPENMM": {
            "format": "AMBER", "topology": "system.prmtop", "coordinates": "system.inpcrd",
        },
        "GROMACS": {
            "format": "GROMACS", "topology": "system.top",
            "topology_hmr": "system_hmr.top", "coordinates": "system.g96",
        },
        "AMBER": {
            "format": "AMBER", "topology": "system.prmtop",
            "topology_hmr": "system_hmr.prmtop", "coordinates": "system.inpcrd",
        },
        "NAMD": {
            "format": "AMBER", "topology": "system.prmtop",
            "topology_hmr": "system_hmr.prmtop", "coordinates": "system.inpcrd",
            "parameters": [],
        },
    }
    return assets, normal_system, runtime_hmr_system


def numerical_amber_gromacs_validation(
    directory: Path,
    reference_variant: str = "native_amber_jac_tip3p",
    cutoff_nm: float = 0.9,
    ewald_error_tolerance: float = 1.0e-5,
) -> dict:
    """Compare native Amber and converted GROMACS systems on one CPU thread."""
    from openmm import Context, Platform, Vec3, VerletIntegrator, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, GromacsTopFile, PME

    amber_top = AmberPrmtopFile(str(directory / "system.prmtop"))
    amber_coordinates = AmberInpcrdFile(str(directory / "system.inpcrd"))
    amber_top.topology.setPeriodicBoxVectors(amber_coordinates.boxVectors)
    _, g96_coordinates_angstrom, g96_box_angstrom = read_g96(directory / "system.g96")
    g96_positions = [Vec3(*(value / 10.0 for value in xyz))
                     for xyz in g96_coordinates_angstrom] * unit.nanometer
    g96_box_vectors = (
        Vec3(g96_box_angstrom[0] / 10.0, 0.0, 0.0),
        Vec3(0.0, g96_box_angstrom[1] / 10.0, 0.0),
        Vec3(0.0, 0.0, g96_box_angstrom[2] / 10.0),
    ) * unit.nanometer
    gromacs_top = GromacsTopFile(
        str(directory / "system.top"), periodicBoxVectors=g96_box_vectors,
        includeDir=str(directory), defines={"FLEXIBLE": "1"},
    )
    options = dict(
        nonbondedMethod=PME, nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=None, rigidWater=False, removeCMMotion=False,
        ewaldErrorTolerance=ewald_error_tolerance,
    )
    systems = (amber_top.createSystem(**options), gromacs_top.createSystem(**options))
    platform = Platform.getPlatformByName("CPU")
    energies = []
    forces = []
    lock_path = Path("/tmp/bench.lock")
    with lock_path.open("a+") as lock_stream:
        fcntl.flock(lock_stream, fcntl.LOCK_EX)
        try:
            for system, positions in zip(systems, (amber_coordinates.positions, g96_positions)):
                integrator = VerletIntegrator(0.001 * unit.picoseconds)
                context = Context(system, integrator, platform, {"Threads": "1"})
                context.setPositions(positions)
                state = context.getState(getEnergy=True, getForces=True)
                energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
                forces.append(state.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoule_per_mole / unit.nanometer
                ))
                del context, integrator
        finally:
            fcntl.flock(lock_stream, fcntl.LOCK_UN)
    numerator = denominator = 0.0
    for expected, observed in zip(forces[0], forces[1]):
        for dimension in range(3):
            delta = float(observed[dimension]) - float(expected[dimension])
            numerator += delta * delta
            denominator += float(expected[dimension]) ** 2
    relative_rms = math.sqrt(numerator / denominator)
    energy_delta = energies[1] - energies[0]
    if abs(energy_delta) > 1.0e-2 or relative_rms > 1.0e-6:
        raise PreparationError(
            "native Amber->GROMACS validation failed: dE=%.9g kJ/mol force RMS=%.9g"
            % (energy_delta, relative_rms)
        )
    return {
        "status": "pass",
        "reference_variant": reference_variant,
        # CPU PME reductions vary at harmless last digits with machine load.
        # These resolutions remain 10x and 100x tighter than the acceptance
        # limits while making independently prepared manifests byte-stable.
        "energy_delta_kj_mol": round(energy_delta, 3),
        "force_relative_rms": round(relative_rms, 8),
        "cutoff_nm": cutoff_nm,
        "ewald_error_tolerance": ewald_error_tolerance,
        "platform": (
            "OpenMM CPU (Threads=1; rounded: deltaE 0.001 kJ/mol, "
            "force RMS 1e-8)"
        ),
    }


def validate_native_dhfr(
    directory: Path,
    neutral: object,
    normal_system: object,
    hmr_system: object,
    hmr_metadata: dict,
    preservation: dict,
) -> dict:
    import parmed as pmd
    from openmm import NonbondedForce, Vec3, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, GromacsTopFile, HBonds, PME

    amber = pmd.load_file(str(directory / "system.prmtop"), xyz=str(directory / "system.inpcrd"))
    amber_hmr = pmd.load_file(str(directory / "system_hmr.prmtop"))
    gromacs = pmd.load_file(str(directory / "system.top"))
    gromacs_hmr = pmd.load_file(str(directory / "system_hmr.top"))
    g96_names, g96_coordinates, g96_box = read_g96(directory / "system.g96")
    openmm_topology = AmberPrmtopFile(str(directory / "system.prmtop"))
    openmm_coordinates = AmberInpcrdFile(str(directory / "system.inpcrd"))
    if openmm_coordinates.boxVectors is None:
        raise PreparationError("OpenMM reload lost the native coordinate box")
    openmm_topology.topology.setPeriodicBoxVectors(openmm_coordinates.boxVectors)
    g96_box_vectors = (
        Vec3(g96_box[0] / 10.0, 0.0, 0.0),
        Vec3(0.0, g96_box[1] / 10.0, 0.0),
        Vec3(0.0, 0.0, g96_box[2] / 10.0),
    ) * unit.nanometer
    gromacs_topology = GromacsTopFile(
        str(directory / "system.top"), periodicBoxVectors=g96_box_vectors,
        includeDir=str(directory),
    )
    gromacs_hmr_topology = GromacsTopFile(
        str(directory / "system_hmr.top"), periodicBoxVectors=g96_box_vectors,
        includeDir=str(directory),
    )
    gromacs_system = gromacs_topology.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds, rigidWater=True, removeCMMotion=False,
    )
    gromacs_hmr_system = gromacs_hmr_topology.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds, rigidWater=True, removeCMMotion=False,
    )

    verify_source_preservation(amber, preservation)
    atom_count = len(neutral.atoms)
    reference_coordinates = neutral.coordinates
    reference_names = [(atom.residue.name, atom.name) for atom in neutral.atoms]
    for label, structure in (
        ("Amber", amber), ("GROMACS", gromacs), ("GROMACS HMR", gromacs_hmr),
    ):
        if len(structure.atoms) != atom_count:
            raise PreparationError("%s reload changed native atom count" % label)
        if [(atom.residue.name, atom.name) for atom in structure.atoms] != reference_names:
            raise PreparationError("%s reload changed native atom/residue order" % label)
        if integer_charge(math.fsum(float(atom.charge) for atom in structure.atoms), label) != 0:
            raise PreparationError("%s reload is not formally neutral" % label)
    if g96_names != reference_names:
        raise PreparationError("G96 reload changed native atom/residue order")
    openmm_names = [
        (atom.residue.name, atom.name) for atom in openmm_topology.topology.atoms()
    ]
    gromacs_openmm_names = [
        (atom.residue.name, atom.name) for atom in gromacs_topology.topology.atoms()
    ]
    gromacs_hmr_openmm_names = [
        (atom.residue.name, atom.name) for atom in gromacs_hmr_topology.topology.atoms()
    ]
    # AmberPrmtopFile applies its documented canonical naming (for example,
    # HID->HIS and terminal H1->H).  The two OpenMM loaders must nevertheless
    # agree index-for-index, while ParmEd/G96 retain the archive names above.
    if openmm_names != gromacs_openmm_names or openmm_names != gromacs_hmr_openmm_names:
        raise PreparationError("OpenMM cross-format topology order/names differ")
    if (openmm_topology.topology.getNumAtoms() != atom_count
            or normal_system.getNumParticles() != atom_count
            or hmr_system.getNumParticles() != atom_count
            or gromacs_system.getNumParticles() != atom_count
            or gromacs_hmr_system.getNumParticles() != atom_count):
        raise PreparationError("OpenMM native asset reload changed particle count")

    if len(amber_hmr.atoms) != atom_count:
        raise PreparationError("static HMR topology changed native atom count")
    for normal_atom, hmr_atom in zip(amber.atoms, amber_hmr.atoms):
        normal_fields = (
            normal_atom.name, normal_atom.residue.name, normal_atom.type,
            normal_atom.atomic_number, float(normal_atom.charge), normal_atom.nb_idx,
        )
        hmr_fields = (
            hmr_atom.name, hmr_atom.residue.name, hmr_atom.type,
            hmr_atom.atomic_number, float(hmr_atom.charge), hmr_atom.nb_idx,
        )
        if normal_fields != hmr_fields:
            raise PreparationError("static HMR changed a non-mass atom parameter")
    if term_connectivity(amber, set()) != term_connectivity(amber_hmr, set()):
        raise PreparationError("static HMR changed bonded connectivity")

    sodium_atoms = [atom for atom in amber.atoms if atom.residue.name == "Na+"]
    if len(sodium_atoms) != preservation["removed_hydrogen_count"] // 2:
        raise PreparationError("reloaded topology has the wrong number of neutralizing Na+ ions")
    for sodium in sodium_atoms:
        if (len(sodium.residue.atoms) != 1 or sodium.name != "Na+"
                or sodium.atomic_number != 11 or abs(float(sodium.charge) - 1.0) > 1.0e-12
                or abs(float(sodium.mass) - preservation["sodium_mass"]) > 1.0e-8
                or abs(float(amber.LJ_radius[sodium.nb_idx - 1])
                       - preservation["sodium_radius"]) > 1.0e-8
                or abs(float(amber.LJ_depth[sodium.nb_idx - 1])
                       - preservation["sodium_epsilon"]) > 1.0e-8
                or sodium.bonds):
            raise PreparationError("reloaded Na+ does not match native JC TIP3P parameters")

    amber_delta = coordinate_delta(reference_coordinates, amber.coordinates)
    gromacs_delta = coordinate_delta(reference_coordinates, g96_coordinates)
    openmm_delta = coordinate_delta(
        reference_coordinates, openmm_coordinates.positions.value_in_unit(unit.angstrom)
    )
    if amber_delta > 1.0e-5 or gromacs_delta > 5.1e-9 or openmm_delta > 1.0e-5:
        raise PreparationError("native coordinate round-trip exceeded file precision")
    box = [float(value) for value in amber.box[:3]]
    if any(abs(box[index] - float(neutral.box[index])) > 1.0e-5 for index in range(3)):
        raise PreparationError("native periodic box changed")
    if any(abs(g96_box[index] - box[index]) > 5.1e-9 for index in range(3)):
        raise PreparationError("G96 periodic box changed")
    for label, candidate_box in (
        ("OpenMM", system_box_angstrom(normal_system)),
        ("OpenMM HMR", system_box_angstrom(hmr_system)),
        ("GROMACS", system_box_angstrom(gromacs_system)),
        ("GROMACS HMR", system_box_angstrom(gromacs_hmr_system)),
    ):
        if any(abs(candidate_box[index] - box[index]) > 1.0e-6 for index in range(3)):
            raise PreparationError("%s System periodic box changed" % label)

    amber_charge = math.fsum(float(atom.charge) for atom in amber.atoms)
    gromacs_charge = math.fsum(float(atom.charge) for atom in gromacs.atoms)
    gromacs_hmr_charge = math.fsum(float(atom.charge) for atom in gromacs_hmr.atoms)
    openmm_charge = system_charge(normal_system)
    gromacs_system_charge = system_charge(gromacs_system)
    gromacs_hmr_system_charge = system_charge(gromacs_hmr_system)
    for charge in (
        amber_charge, gromacs_charge, gromacs_hmr_charge, openmm_charge,
        gromacs_system_charge, gromacs_hmr_system_charge,
    ):
        if integer_charge(charge, "native reloaded asset") != 0:
            raise PreparationError("native reloaded asset is not neutral")
    reference_signature = system_topology_signature(normal_system)
    # ParmEd intentionally omits zero-force Amber torsion records from the
    # GROMACS topology.  They contribute no bonded energy; their exclusion
    # bookkeeping is independently protected by the exact exception-pair set.
    for label, signature in (
        ("OpenMM HMR", system_topology_signature(hmr_system)),
        ("GROMACS", system_topology_signature(gromacs_system)),
        ("GROMACS HMR", system_topology_signature(gromacs_hmr_system)),
    ):
        if signature["bonded_force_term_counts"] != reference_signature["bonded_force_term_counts"]:
            raise PreparationError("%s bonded force term counts changed" % label)
        if signature["exception_pairs"] != reference_signature["exception_pairs"]:
            raise PreparationError("%s exception/1-4 pair set changed" % label)
        if set(signature["constraints"]) != set(reference_signature["constraints"]):
            raise PreparationError("%s constraint pair set changed" % label)
        if signature["virtual_sites"] != reference_signature["virtual_sites"]:
            raise PreparationError("%s virtual-site set changed" % label)
        if any(abs(signature["constraints"][pair] - distance) > 1.0e-8
               for pair, distance in reference_signature["constraints"].items()):
            raise PreparationError("%s constraint distances changed" % label)
        if label == "OpenMM HMR" and (
                signature["zero_energy_torsions"]
                != reference_signature["zero_energy_torsions"]):
            raise PreparationError("OpenMM HMR zero-energy torsion count changed")
        if label.startswith("GROMACS") and signature["zero_energy_torsions"] != 0:
            raise PreparationError("GROMACS export retained an unexpected zero-energy torsion")
    normal_mass = mass_structure(amber)
    amber_hmr_mass = mass_structure(amber_hmr)
    gromacs_mass = mass_structure(gromacs)
    gromacs_hmr_mass = mass_structure(gromacs_hmr)
    openmm_mass = mass_system(normal_system)
    openmm_hmr_mass = mass_system(hmr_system)
    if (abs(amber_hmr_mass - normal_mass) > 1.0e-6
            or abs(gromacs_mass - normal_mass) > 1.0e-4
            or abs(gromacs_hmr_mass - normal_mass) > 1.0e-4):
        raise PreparationError("native static HMR/GROMACS export changed total mass")
    if abs(openmm_hmr_mass - openmm_mass) > 1.0e-6:
        raise PreparationError("native OpenMM HMR export changed total mass")
    normal_particle_masses = [float(atom.mass) for atom in amber.atoms]
    static_hmr_particle_masses = [float(atom.mass) for atom in amber_hmr.atoms]
    runtime_hmr_particle_masses = [
        hmr_system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(atom_count)
    ]
    if any(abs(expected - observed) > 1.0e-6 for expected, observed in zip(
            static_hmr_particle_masses, runtime_hmr_particle_masses)):
        raise PreparationError("static and OpenMM runtime HMR particle masses differ")
    if any(abs(expected - float(observed.mass)) > 1.0e-6 for expected, observed in zip(
            static_hmr_particle_masses, gromacs_hmr.atoms)):
        raise PreparationError("static Amber/GROMACS HMR particle masses differ")
    changed_particles = {
        index for index, (before, after) in enumerate(zip(
            normal_particle_masses, static_hmr_particle_masses
        )) if abs(before - after) > 1.0e-8
    }
    if len(changed_particles) != hmr_metadata["particles_with_changed_mass"]:
        raise PreparationError("HMR changed-particle metadata is inconsistent")
    changed_water_particles = sum(
        amber.atoms[index].residue.name.upper() in WATER_NAMES for index in changed_particles
    )
    if changed_water_particles != hmr_metadata["water_particles_with_changed_mass"]:
        raise PreparationError("HMR unexpectedly changed a rigid-water particle")

    impropers = sum(bool(getattr(item, "improper", False)) for item in amber.dihedrals)
    nonbonded = [force for force in normal_system.getForces() if isinstance(force, NonbondedForce)]
    canonical = {
        "atom_count": atom_count,
        "raw_charge_e": amber_charge,
        "formal_charge_e": integer_charge(amber_charge, "native canonical asset"),
        "box_angstrom": box,
        "bonded_term_counts": {
            "bonds": len(amber.bonds), "angles": len(amber.angles),
            "proper_dihedrals": len(amber.dihedrals) - impropers,
            "impropers": len(amber.impropers) + impropers, "cmaps": len(amber.cmaps),
        },
        "exception_count": nonbonded[0].getNumExceptions(),
        "constraint_count": normal_system.getNumConstraints(),
    }
    amber_record = asset_record(
        atom_count, amber_charge, box, amber_delta, normal_mass, amber_hmr_mass
    )
    assets = {
        "GENESIS": dict(amber_record, hmr_total_mass_dalton=None),
        "OPENMM": asset_record(
            atom_count, openmm_charge, box, openmm_delta, openmm_mass, openmm_hmr_mass,
        ),
        "GROMACS": asset_record(
            atom_count, gromacs_charge, g96_box, gromacs_delta,
            gromacs_mass, gromacs_hmr_mass,
        ),
        "AMBER": dict(amber_record),
        "NAMD": dict(amber_record),
    }
    numerical = numerical_amber_gromacs_validation(directory)
    return {"status": "pass", "canonical": canonical, "assets": assets,
            "representative_numerical": numerical}


def amber_flag_differences(
    source_path: Path,
    candidate_path: Path,
    allowed: set[str],
) -> set[str]:
    """Require two prmtops to differ only in explicitly allowed arrays."""
    import parmed as pmd

    source = pmd.load_file(str(source_path))
    candidate = pmd.load_file(str(candidate_path))
    flags = set(source.parm_data) | set(candidate.parm_data)
    differences = {
        flag for flag in flags
        if flag not in source.parm_data or flag not in candidate.parm_data
        or list(source.parm_data[flag]) != list(candidate.parm_data[flag])
    }
    unexpected = differences - allowed
    if unexpected:
        raise PreparationError(
            "%s changed native Amber arrays outside %s: %s" % (
                candidate_path.name, ", ".join(sorted(allowed)),
                ", ".join(sorted(unexpected)),
            )
        )
    return differences


def reconcile_native_amber_box(topology: Path, coordinates: Path) -> dict:
    """Patch only a stale prmtop BOX_DIMENSIONS field from its restart."""
    import parmed as pmd

    topology_only = pmd.load_file(str(topology))
    with_coordinates = pmd.load_file(str(topology), xyz=str(coordinates))
    if topology_only.box is None or with_coordinates.box is None:
        raise PreparationError("native Amber explicit-solvent asset has no periodic box")
    source_box = [float(value) for value in topology_only.box]
    coordinate_box = [float(value) for value in with_coordinates.box]
    changed = any(abs(left - right) > 1.0e-10 for left, right in zip(
        source_box, coordinate_box,
    ))
    if not changed:
        return {
            "patched": False,
            "source_box_angstrom_degrees": source_box,
            "coordinate_box_angstrom_degrees": coordinate_box,
        }
    if any(abs(angle - 90.0) > 1.0e-10 for angle in coordinate_box[3:6]):
        raise PreparationError(
            "traditional Amber BOX_DIMENSIONS cannot encode this restart box"
        )

    lines = topology.read_text(encoding="ascii").splitlines(keepends=True)
    try:
        flag = next(
            index for index, line in enumerate(lines)
            if line.strip() == "%FLAG BOX_DIMENSIONS"
        )
    except StopIteration as error:
        raise PreparationError("periodic Amber topology has no BOX_DIMENSIONS") from error
    if flag + 2 >= len(lines) or lines[flag + 1].strip().upper() != "%FORMAT(5E16.8)":
        raise PreparationError("unsupported Amber BOX_DIMENSIONS format")
    end = next(
        (index for index in range(flag + 2, len(lines)) if lines[index].startswith("%FLAG")),
        len(lines),
    )
    values = [coordinate_box[4], *coordinate_box[:3]]
    serialized_values = [float("%.8E" % value) for value in values]
    lines[flag + 2:end] = ["".join("%16.8E" % value for value in values) + "\n"]
    topology.write_text("".join(lines), encoding="ascii")

    reloaded = pmd.load_file(str(topology))
    if reloaded.box is None or any(
            abs(float(observed) - expected) > 1.0e-12
            for observed, expected in zip(
                reloaded.box,
                [serialized_values[1], serialized_values[2], serialized_values[3],
                 serialized_values[0], serialized_values[0], serialized_values[0]],
            )):
        raise PreparationError("patched Amber BOX_DIMENSIONS failed round-trip validation")
    return {
        "patched": True,
        "scope": "BOX_DIMENSIONS only",
        "source_box_angstrom_degrees": source_box,
        "coordinate_box_angstrom_degrees": coordinate_box,
    }


def canonicalize_native_amber_multiterm_dihedrals(topology: Path) -> dict:
    """Expand Amber's negative-periodicity continuation representation.

    Amber uses a negative ``DIHEDRAL_PERIODICITY`` value to say that the next
    parameter-table entry is another Fourier term for the same atom quartet.
    pmemd expands those entries during setup and then takes the absolute value
    of every periodicity.  OpenMM's direct prmtop reader and ParmEd's topology
    converter do not expand that representation.  Make the implicit terms
    explicit, preserving exactly one 1-4 interaction per original quartet by
    setting Amber's negative-third-atom suppression flag on nonterminal terms.
    """
    from parmed.amber import AmberFormat, AmberParm
    from parmed.constants import PrmtopPointers

    amber = AmberParm(str(topology))
    periodicities = list(amber.parm_data["DIHEDRAL_PERIODICITY"])
    if any(not math.isfinite(float(value)) for value in periodicities):
        raise PreparationError("native Amber topology has a nonfinite periodicity")
    if any(abs(float(value) - round(float(value))) > 1.0e-8
           for value in periodicities):
        raise PreparationError("native Amber topology has a nonintegral periodicity")
    if any(float(value) == 0.0 for value in periodicities):
        raise PreparationError(
            "zero-periodicity native Amber torsions require a separate exact normalization"
        )
    negative_types = [
        index for index, value in enumerate(periodicities, 1) if float(value) < 0.0
    ]
    if not negative_types:
        return {
            "expanded": False,
            "changed_flags": set(),
            "negative_parameter_type_count": 0,
            "original_term_counts": {},
            "expanded_term_counts": {},
        }

    changed_flags = {"DIHEDRAL_PERIODICITY", "POINTERS"}
    original_counts = {}
    expanded_counts = {}
    arrays = (
        ("DIHEDRALS_INC_HYDROGEN", PrmtopPointers.NPHIH, None),
        ("DIHEDRALS_WITHOUT_HYDROGEN", PrmtopPointers.NPHIA,
         PrmtopPointers.MPHIA),
    )
    for flag, primary_pointer, matching_pointer in arrays:
        raw = list(amber.parm_data[flag])
        if len(raw) % 5:
            raise PreparationError("malformed native Amber %s array" % flag)
        expanded = []
        for offset in range(0, len(raw), 5):
            atom1, atom2, atom3, atom4, type_index = raw[offset:offset + 5]
            sequence = []
            current = type_index
            while True:
                if current <= 0 or current > len(periodicities):
                    raise PreparationError(
                        "native Amber multiterm dihedral runs outside its parameter table"
                    )
                sequence.append(current)
                if periodicities[current - 1] > 0.0:
                    break
                current += 1
            if atom3 == 0 and len(sequence) > 1:
                raise PreparationError(
                    "cannot encode 1-4 suppression for a multiterm dihedral at atom zero"
                )
            for sequence_index, current in enumerate(sequence):
                candidate_atom3 = atom3
                if atom3 > 0 and sequence_index < len(sequence) - 1:
                    candidate_atom3 = -atom3
                expanded.extend((
                    atom1, atom2, candidate_atom3, atom4, current,
                ))
        original_counts[flag] = len(raw) // 5
        expanded_counts[flag] = len(expanded) // 5
        amber.parm_data[flag] = expanded
        amber.parm_data["POINTERS"][primary_pointer] = len(expanded) // 5
        if matching_pointer is not None:
            amber.parm_data["POINTERS"][matching_pointer] = len(expanded) // 5
        if expanded != raw:
            changed_flags.add(flag)

    amber.parm_data["DIHEDRAL_PERIODICITY"] = [
        abs(float(value)) for value in periodicities
    ]
    # AmberParm.write_parm calls remake_parm(), which would reconstruct the
    # original unexpanded Structure.  Invoke the raw AmberFormat writer after
    # changing the authoritative arrays instead.
    AmberFormat.write_parm(amber, str(topology))
    normalize_amber_header(topology)

    reloaded = AmberParm(str(topology))
    if any(float(value) <= 0.0
           for value in reloaded.parm_data["DIHEDRAL_PERIODICITY"]):
        raise PreparationError("canonical Amber topology retained a nonpositive periodicity")
    for flag, primary_pointer, matching_pointer in arrays:
        expected = expanded_counts[flag]
        if (len(reloaded.parm_data[flag]) != expected * 5
                or reloaded.ptr(primary_pointer.name) != expected
                or (matching_pointer is not None
                    and reloaded.ptr(matching_pointer.name) != expected)):
            raise PreparationError("canonical Amber multiterm expansion did not round-trip")
    return {
        "expanded": True,
        "changed_flags": changed_flags,
        "negative_parameter_type_count": len(negative_types),
        "original_term_counts": original_counts,
        "expanded_term_counts": expanded_counts,
    }


def native_amber_source_transformations(
    box_reconciliation: dict,
    dihedral_normalization: dict,
) -> dict:
    """Serialize every sanctioned source-prmtop canonicalization."""
    box = {
        "patched": bool(box_reconciliation["patched"]),
        "scope": box_reconciliation.get("scope"),
        "source_box_angstrom_degrees": list(
            box_reconciliation["source_box_angstrom_degrees"]
        ),
        "coordinate_box_angstrom_degrees": list(
            box_reconciliation["coordinate_box_angstrom_degrees"]
        ),
    }
    dihedrals = {
        "expanded": bool(dihedral_normalization["expanded"]),
        "changed_flags": sorted(dihedral_normalization["changed_flags"]),
        "negative_parameter_type_count": int(
            dihedral_normalization["negative_parameter_type_count"]
        ),
        "original_term_counts": dict(sorted(
            dihedral_normalization["original_term_counts"].items()
        )),
        "expanded_term_counts": dict(sorted(
            dihedral_normalization["expanded_term_counts"].items()
        )),
    }
    return {
        "amber_box_dimensions": box,
        "amber_multiterm_dihedrals": dihedrals,
    }


def native_amber_none_required(structure: object) -> dict:
    """Describe a formally neutral native AMBER system without mutation."""
    raw_charge = math.fsum(float(atom.charge) for atom in structure.atoms)
    if integer_charge(raw_charge, "native Amber source") != 0:
        raise PreparationError("native Amber source unexpectedly requires counterions")
    water_count = 0
    existing_ions: dict[str, int] = {}
    ion_atomic_numbers = {11, 12, 17, 19, 20}
    for residue in structure.residues:
        atoms = list(residue.atoms)
        elements = sorted(atom.atomic_number for atom in atoms)
        if residue.name.upper() in WATER_NAMES and elements == [1, 1, 8]:
            water_count += 1
        elif len(atoms) == 1 and atoms[0].atomic_number in ion_atomic_numbers:
            existing_ions[residue.name] = existing_ions.get(residue.name, 0) + 1
    return {
        "method": "none_required",
        "seed": None,
        "tolerance_e": FORMAL_CHARGE_TOLERANCE_E,
        "min_ion_separation_nm": None,
        "min_ion_solute_distance_nm": None,
        "source_charge_e": raw_charge,
        "source_formal_charge_e": 0,
        "pre_neutralization_charge_e": raw_charge,
        "pre_neutralization_formal_charge_e": 0,
        "post_neutralization_charge_e": raw_charge,
        "post_neutralization_formal_charge_e": 0,
        "neutralization_ion": None,
        "neutralization_ion_charge_e": None,
        "neutralization_ion_count": 0,
        "replaced_water_count": 0,
        "existing_ions": dict(sorted(existing_ions.items())),
        "water_count_before": water_count,
        "water_count_after": water_count,
        "selected_water_residue_indices": [],
        "selected_water_oxygen_positions_nm": [],
    }


def write_native_amber_gromacs_topology(structure: object, path: Path) -> None:
    """Write an inline topology with only physically active 1-4 scales."""
    import parmed as pmd
    from parmed.topologyobjects import DihedralTypeList

    converted = pmd.gromacs.GromacsTopologyFile.from_structure(structure)
    active_scee = set()
    active_scnb = set()
    for dihedral in structure.dihedrals:
        if dihedral.ignore_end or dihedral.type is None:
            continue
        types = (
            dihedral.type if isinstance(dihedral.type, DihedralTypeList)
            else (dihedral.type,)
        )
        for term in types:
            if (abs(float(dihedral.atom1.charge) * float(dihedral.atom4.charge)) > 1.0e-16):
                if term.scee is None or float(term.scee) <= 0.0:
                    raise PreparationError("active Amber 1-4 electrostatics has invalid SCEE")
                active_scee.add(float(term.scee))
            if (float(dihedral.atom1.epsilon) > 0.0
                    and float(dihedral.atom4.epsilon) > 0.0):
                if term.scnb is None or float(term.scnb) <= 0.0:
                    raise PreparationError("active Amber 1-4 Lennard-Jones has invalid SCNB")
                active_scnb.add(float(term.scnb))

    for label, values in (("SCEE", active_scee), ("SCNB", active_scnb)):
        if len({round(value, 8) for value in values}) > 1:
            raise PreparationError(
                "mixed active Amber %s factors are not representable in GROMACS: %s" % (
                    label, ", ".join("%.12g" % value for value in sorted(values)),
                )
            )
    if active_scee:
        converted.defaults.fudgeQQ = 1.0 / next(iter(active_scee))
    if active_scnb:
        converted.defaults.fudgeLJ = 1.0 / next(iter(active_scnb))
    converted.write(str(path), combine=None, parameters="inline")
    normalize_gromacs_header(path)
    rewrite_gromacs_nonbonded_precision(converted, path)


def rewrite_gromacs_nonbonded_precision(topology: object, path: Path) -> None:
    """Restore source precision lost by ParmEd's GROMACS text formatter.

    ParmEd writes charges with eight decimal places and atom-type Lennard-Jones
    values with eight significant digits.  That is enough for ordinary use but
    can move a large-system PME single point beyond this benchmark's 0.01
    kJ/mol equivalence threshold.  Rewrite only the defaults, atom types, and
    per-molecule atom charge/mass fields from the in-memory converted topology.
    """
    from openmm import unit
    from parmed.parameters import ParameterSet

    parameters = ParameterSet.from_structure(
        topology, allow_unequal_duplicates=True,
    )
    molecules = [molecule for molecule, _ in topology.split()]
    print_bond_types = any(
        atom_type._bond_type is not None
        for atom_type in parameters.atom_types.values()
    )
    print_atomic_numbers = all(
        atom_type.atomic_number != -1
        for atom_type in parameters.atom_types.values()
    )
    energy_conversion = unit.kilocalories.conversion_factor_to(unit.kilojoules)

    def precise(value: float) -> str:
        return format(float(value), ".17e")

    section = None
    molecule_index = -1
    atom_rows = Counter()
    atom_type_rows = set()
    rewritten = []
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip()
            rewritten.append(line)
            continue
        body, separator, comment = line.partition(";")
        fields = body.split()
        is_data = bool(fields and not fields[0].startswith("#"))
        if is_data and section == "moleculetype":
            molecule_index += 1
            if molecule_index >= len(molecules):
                raise PreparationError("GROMACS writer emitted an extra molecule type")
        if is_data and section == "defaults":
            if len(fields) < 5:
                raise PreparationError("GROMACS defaults row lost its scaling factors")
            fields[3] = precise(topology.defaults.fudgeLJ)
            fields[4] = precise(topology.defaults.fudgeQQ)
            body = " ".join(fields)
        elif is_data and section == "atomtypes":
            try:
                atom_type = parameters.atom_types[fields[0]]
            except KeyError as error:
                raise PreparationError(
                    "GROMACS atom type %s was not present in the converted topology"
                    % fields[0]
                ) from error
            offset = 1 + int(print_bond_types) + int(print_atomic_numbers)
            if len(fields) != offset + 5:
                raise PreparationError("unexpected generated GROMACS atom-type layout")
            fields[offset] = precise(atom_type.mass)
            fields[offset + 1] = precise(atom_type.charge)
            fields[offset + 3] = precise(atom_type.sigma / 10.0)
            fields[offset + 4] = precise(atom_type.epsilon * energy_conversion)
            body = " ".join(fields)
            atom_type_rows.add(fields[0])
        elif is_data and section == "atoms":
            if molecule_index < 0:
                raise PreparationError("GROMACS atoms row precedes its molecule type")
            try:
                atom = molecules[molecule_index].atoms[int(fields[0]) - 1]
            except (IndexError, ValueError) as error:
                raise PreparationError("invalid generated GROMACS atom index") from error
            if len(fields) < 8:
                raise PreparationError("generated GROMACS atom row lost charge or mass")
            fields[6] = precise(atom.charge)
            fields[7] = precise(atom.mass)
            body = " ".join(fields)
            atom_rows[molecule_index] += 1
        rewritten.append(body + (separator + comment if separator else ""))

    if molecule_index + 1 != len(molecules):
        raise PreparationError("GROMACS writer omitted a molecule type")
    if any(atom_rows[index] != len(molecule.atoms)
           for index, molecule in enumerate(molecules)):
        raise PreparationError("GROMACS high-precision rewrite missed an atom row")
    if atom_type_rows != set(parameters.atom_types):
        raise PreparationError("GROMACS high-precision rewrite missed an atom type")
    path.write_text("\n".join(rewritten) + "\n", encoding="ascii")


def compare_system_signature(reference: dict, candidate: dict, label: str) -> None:
    """Require a converted or HMR System to preserve active physical terms."""
    if candidate["bonded_force_term_counts"] != reference["bonded_force_term_counts"]:
        raise PreparationError("%s bonded force term counts changed" % label)
    if candidate["exception_pairs"] != reference["exception_pairs"]:
        raise PreparationError("%s exception/1-4 pair set changed" % label)
    if set(candidate["constraints"]) != set(reference["constraints"]):
        raise PreparationError("%s constraint pair set changed" % label)
    if candidate["virtual_sites"] != reference["virtual_sites"]:
        raise PreparationError("%s virtual-site set changed" % label)
    if any(
            abs(candidate["constraints"][pair] - distance) > 1.0e-8
            for pair, distance in reference["constraints"].items()):
        raise PreparationError("%s constraint distances changed" % label)


def export_native_amber_archive(
    directory: Path,
    normal: object,
    hmr: object,
    include_gromacs: bool,
) -> tuple[dict, object, object]:
    """Export direct native-AMBER assets and proven GROMACS overlays."""
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME

    hmr.save(str(directory / "system_hmr.prmtop"), overwrite=False)
    normalize_amber_header(directory / "system_hmr.prmtop")

    amber_topology = AmberPrmtopFile(str(directory / "system.prmtop"))
    coordinates = AmberInpcrdFile(str(directory / "system.inpcrd"))
    if coordinates.boxVectors is None:
        raise PreparationError("written native Amber restart has no periodic box")
    amber_topology.topology.setPeriodicBoxVectors(coordinates.boxVectors)
    options = dict(
        nonbondedMethod=PME,
        nonbondedCutoff=0.9 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=False,
    )
    normal_system = amber_topology.createSystem(**options)
    runtime_hmr_system = amber_topology.createSystem(
        **options, hydrogenMass=HMR_HYDROGEN_MASS_DA * unit.dalton,
    )

    assets = {
        "GENESIS": {
            "format": "AMBER", "topology": "system.prmtop",
            "coordinates": "system.inpcrd", "topology_definitions": [], "parameters": [],
        },
        "OPENMM": {
            "format": "AMBER", "topology": "system.prmtop",
            "coordinates": "system.inpcrd",
        },
        "AMBER": {
            "format": "AMBER", "topology": "system.prmtop",
            "topology_hmr": "system_hmr.prmtop", "coordinates": "system.inpcrd",
        },
        "NAMD": {
            "format": "AMBER", "topology": "system.prmtop",
            "topology_hmr": "system_hmr.prmtop", "coordinates": "system.inpcrd",
            "parameters": [],
        },
    }
    if include_gromacs:
        write_native_amber_gromacs_topology(normal, directory / "system.top")
        write_native_amber_gromacs_topology(hmr, directory / "system_hmr.top")
        write_g96(directory / "system.g96", normal)
        assets["GROMACS"] = {
            "format": "GROMACS", "topology": "system.top",
            "topology_hmr": "system_hmr.top", "coordinates": "system.g96",
        }
    return assets, normal_system, runtime_hmr_system


def validate_native_amber_archive(
    directory: Path,
    source_topology: Path,
    source_coordinates: Path,
    normal: object,
    normal_system: object,
    runtime_hmr_system: object,
    hmr_metadata: dict,
    assets: dict,
    box_reconciliation: dict,
    dihedral_normalization: dict,
    reference_variant: str,
) -> dict:
    """Reload and validate every advertised native-AMBER engine asset."""
    import parmed as pmd
    from openmm import NonbondedForce, Vec3, unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, GromacsTopFile, HBonds, PME

    source = pmd.load_file(str(source_topology), xyz=str(source_coordinates))
    amber = pmd.load_file(
        str(directory / "system.prmtop"), xyz=str(directory / "system.inpcrd")
    )
    amber_hmr = pmd.load_file(str(directory / "system_hmr.prmtop"))
    normal_changes = amber_flag_differences(
        source_topology, directory / "system.prmtop", {
            "BOX_DIMENSIONS", "POINTERS", "DIHEDRAL_PERIODICITY",
            "DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN",
        },
    )
    expected_normal_changes = set(dihedral_normalization["changed_flags"])
    if box_reconciliation["patched"]:
        expected_normal_changes.add("BOX_DIMENSIONS")
    if normal_changes != expected_normal_changes:
        raise PreparationError("normal Amber overlay did not preserve the source topology")
    hmr_changes = amber_flag_differences(
        directory / "system.prmtop", directory / "system_hmr.prmtop",
        {"MASS", "ATOMIC_NUMBER"},
    )
    if "MASS" not in hmr_changes:
        raise PreparationError("static HMR topology did not change the MASS array")

    atom_count = len(normal.atoms)
    if len(source.atoms) != atom_count or len(amber.atoms) != atom_count:
        raise PreparationError("normal Amber overlay changed atom count")
    reference_names = [(atom.residue.name, atom.name) for atom in normal.atoms]
    for source_atom, candidate in zip(source.atoms, amber.atoms):
        source_fields = (
            source_atom.name, source_atom.residue.name, source_atom.type,
            source_atom.atomic_number, source_atom.nb_idx,
            float(source_atom.charge), float(source_atom.mass),
        )
        candidate_fields = (
            candidate.name, candidate.residue.name, candidate.type,
            candidate.atomic_number, candidate.nb_idx,
            float(candidate.charge), float(candidate.mass),
        )
        if source_fields != candidate_fields:
            raise PreparationError("normal Amber overlay changed an atom parameter")
    if coordinate_delta(source.coordinates, amber.coordinates) > 1.0e-12:
        raise PreparationError("normal Amber overlay changed restart coordinates")
    if len(amber_hmr.atoms) != atom_count:
        raise PreparationError("static HMR topology changed atom count")
    for normal_atom, hmr_atom in zip(amber.atoms, amber_hmr.atoms):
        normal_fields = (
            normal_atom.name, normal_atom.residue.name, normal_atom.type,
            normal_atom.atomic_number, normal_atom.nb_idx, float(normal_atom.charge),
        )
        hmr_fields = (
            hmr_atom.name, hmr_atom.residue.name, hmr_atom.type,
            hmr_atom.atomic_number, hmr_atom.nb_idx, float(hmr_atom.charge),
        )
        if normal_fields != hmr_fields:
            raise PreparationError("static HMR changed a non-mass atom parameter")

    openmm_topology = AmberPrmtopFile(str(directory / "system.prmtop"))
    openmm_coordinates = AmberInpcrdFile(str(directory / "system.inpcrd"))
    if openmm_coordinates.boxVectors is None:
        raise PreparationError("OpenMM reload lost the native Amber restart box")
    openmm_topology.topology.setPeriodicBoxVectors(openmm_coordinates.boxVectors)
    if (openmm_topology.topology.getNumAtoms() != atom_count
            or normal_system.getNumParticles() != atom_count
            or runtime_hmr_system.getNumParticles() != atom_count):
        raise PreparationError("OpenMM reload changed native Amber particle count")
    openmm_delta = coordinate_delta(
        amber.coordinates,
        openmm_coordinates.positions.value_in_unit(unit.angstrom),
    )
    if openmm_delta > 1.0e-5:
        raise PreparationError("OpenMM Amber coordinate reload exceeded restart precision")

    box = [float(value) for value in amber.box[:3]]
    for label, candidate_box in (
        ("OpenMM", system_box_angstrom(normal_system)),
        ("OpenMM HMR", system_box_angstrom(runtime_hmr_system)),
    ):
        if any(abs(left - right) > 1.0e-6 for left, right in zip(candidate_box, box)):
            raise PreparationError("%s System periodic box changed" % label)
    amber_charge = math.fsum(float(atom.charge) for atom in amber.atoms)
    openmm_charge = system_charge(normal_system)
    for label, charge in (("Amber", amber_charge), ("OpenMM", openmm_charge)):
        if integer_charge(charge, label) != 0:
            raise PreparationError("%s native topology is not formally neutral" % label)

    reference_signature = system_topology_signature(normal_system)
    runtime_hmr_signature = system_topology_signature(runtime_hmr_system)
    compare_system_signature(reference_signature, runtime_hmr_signature, "OpenMM HMR")
    if runtime_hmr_signature["zero_energy_torsions"] != reference_signature["zero_energy_torsions"]:
        raise PreparationError("OpenMM HMR zero-energy torsion count changed")

    normal_mass = mass_structure(amber)
    static_hmr_mass = mass_structure(amber_hmr)
    openmm_mass = mass_system(normal_system)
    openmm_hmr_mass = mass_system(runtime_hmr_system)
    if abs(static_hmr_mass - normal_mass) > 1.0e-6:
        raise PreparationError("static HMR changed total native Amber mass")
    if abs(openmm_hmr_mass - openmm_mass) > 1.0e-6:
        raise PreparationError("OpenMM runtime HMR changed total mass")
    static_hmr_masses = [float(atom.mass) for atom in amber_hmr.atoms]
    runtime_hmr_masses = [
        runtime_hmr_system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(atom_count)
    ]
    if any(abs(left - right) > 1.0e-6 for left, right in zip(
            static_hmr_masses, runtime_hmr_masses)):
        raise PreparationError("static and OpenMM runtime HMR masses differ")
    changed_particles = {
        index for index, (before, after) in enumerate(zip(
            (float(atom.mass) for atom in amber.atoms), static_hmr_masses,
        )) if abs(before - after) > 1.0e-8
    }
    if len(changed_particles) != hmr_metadata["particles_with_changed_mass"]:
        raise PreparationError("HMR changed-particle metadata is inconsistent")
    if any(amber.atoms[index].residue.name.upper() in WATER_NAMES
           for index in changed_particles):
        raise PreparationError("HMR changed a rigid-water particle")

    amber_record = asset_record(
        atom_count, amber_charge, box, 0.0, normal_mass, static_hmr_mass,
    )
    validation_assets = {
        "GENESIS": dict(amber_record, hmr_total_mass_dalton=None),
        "OPENMM": asset_record(
            atom_count, openmm_charge, box, openmm_delta, openmm_mass, openmm_hmr_mass,
        ),
        "AMBER": dict(amber_record),
        "NAMD": dict(amber_record),
    }

    numerical = {
        "status": "not_run", "reference_variant": None,
        "energy_delta_kj_mol": None, "force_relative_rms": None,
        "cutoff_nm": None, "ewald_error_tolerance": None, "platform": None,
    }
    if "GROMACS" in assets:
        gromacs = pmd.load_file(str(directory / "system.top"))
        gromacs_hmr = pmd.load_file(str(directory / "system_hmr.top"))
        g96_names, g96_coordinates, g96_box = read_g96(directory / "system.g96")
        if (len(gromacs.atoms) != atom_count or len(gromacs_hmr.atoms) != atom_count
                or [(atom.residue.name, atom.name) for atom in gromacs.atoms] != reference_names
                or g96_names != reference_names):
            raise PreparationError("GROMACS conversion changed native atom order")
        gromacs_delta = coordinate_delta(amber.coordinates, g96_coordinates)
        if gromacs_delta > 5.1e-9:
            raise PreparationError("G96 coordinate round-trip exceeded fixed precision")
        if any(abs(left - right) > 5.1e-9 for left, right in zip(g96_box, box)):
            raise PreparationError("G96 periodic box changed")
        gromacs_charge = math.fsum(float(atom.charge) for atom in gromacs.atoms)
        if integer_charge(gromacs_charge, "GROMACS") != 0:
            raise PreparationError("GROMACS conversion is not formally neutral")
        gromacs_mass = mass_structure(gromacs)
        gromacs_hmr_mass = mass_structure(gromacs_hmr)
        if (abs(gromacs_mass - normal_mass) > 1.0e-4
                or abs(gromacs_hmr_mass - normal_mass) > 1.0e-4):
            raise PreparationError("GROMACS conversion changed total mass")
        if any(abs(float(atom.mass) - expected) > 1.0e-6 for atom, expected in zip(
                gromacs_hmr.atoms, static_hmr_masses)):
            raise PreparationError("GROMACS and Amber static-HMR masses differ")

        g96_box_vectors = (
            Vec3(g96_box[0] / 10.0, 0.0, 0.0),
            Vec3(0.0, g96_box[1] / 10.0, 0.0),
            Vec3(0.0, 0.0, g96_box[2] / 10.0),
        ) * unit.nanometer
        gromacs_options = dict(
            periodicBoxVectors=g96_box_vectors, includeDir=str(directory),
        )
        gromacs_topology = GromacsTopFile(str(directory / "system.top"), **gromacs_options)
        gromacs_hmr_topology = GromacsTopFile(
            str(directory / "system_hmr.top"), **gromacs_options
        )
        create_options = dict(
            nonbondedMethod=PME, nonbondedCutoff=0.9 * unit.nanometer,
            constraints=HBonds, rigidWater=True, removeCMMotion=False,
        )
        gromacs_system = gromacs_topology.createSystem(**create_options)
        gromacs_hmr_system = gromacs_hmr_topology.createSystem(**create_options)
        if (gromacs_system.getNumParticles() != atom_count
                or gromacs_hmr_system.getNumParticles() != atom_count):
            raise PreparationError("GROMACS OpenMM reload changed particle count")
        gromacs_signature = system_topology_signature(gromacs_system)
        gromacs_hmr_signature = system_topology_signature(gromacs_hmr_system)
        compare_system_signature(reference_signature, gromacs_signature, "GROMACS")
        compare_system_signature(reference_signature, gromacs_hmr_signature, "GROMACS HMR")
        if (gromacs_signature["zero_energy_torsions"] != 0
                or gromacs_hmr_signature["zero_energy_torsions"] != 0):
            raise PreparationError("GROMACS retained an unexpected zero-energy torsion")
        for label, candidate_box in (
            ("GROMACS", system_box_angstrom(gromacs_system)),
            ("GROMACS HMR", system_box_angstrom(gromacs_hmr_system)),
        ):
            if any(abs(left - right) > 1.0e-6 for left, right in zip(candidate_box, box)):
                raise PreparationError("%s System periodic box changed" % label)
        if integer_charge(system_charge(gromacs_system), "GROMACS System") != 0:
            raise PreparationError("GROMACS OpenMM System is not formally neutral")
        if (abs(mass_system(gromacs_system) - normal_mass) > 1.0e-4
                or abs(mass_system(gromacs_hmr_system) - normal_mass) > 1.0e-4):
            raise PreparationError("GROMACS OpenMM System changed total mass")
        validation_assets["GROMACS"] = asset_record(
            # The G96 box has already been checked against the authoritative
            # restart above.  Record that canonical box rather than the
            # binary-float result of converting the serialized nm values back
            # to Angstrom (for example, 8.326316990 * 10 can become
            # 83.26316990000001).  Prepared manifests intentionally use one
            # exact box triple for every engine.
            atom_count, gromacs_charge, box, gromacs_delta,
            gromacs_mass, gromacs_hmr_mass,
        )
        numerical = numerical_amber_gromacs_validation(
            directory, reference_variant=reference_variant,
            cutoff_nm=0.9, ewald_error_tolerance=1.0e-5,
        )

    impropers = sum(bool(getattr(item, "improper", False)) for item in amber.dihedrals)
    nonbonded = [
        force for force in normal_system.getForces() if isinstance(force, NonbondedForce)
    ]
    if len(nonbonded) != 1:
        raise PreparationError("native Amber System does not have one NonbondedForce")
    canonical = {
        "atom_count": atom_count,
        "raw_charge_e": amber_charge,
        "formal_charge_e": 0,
        "box_angstrom": box,
        "bonded_term_counts": {
            "bonds": len(amber.bonds), "angles": len(amber.angles),
            "proper_dihedrals": len(amber.dihedrals) - impropers,
            "impropers": len(amber.impropers) + impropers,
            "cmaps": len(amber.cmaps),
        },
        "exception_count": nonbonded[0].getNumExceptions(),
        "constraint_count": normal_system.getNumConstraints(),
    }
    if set(validation_assets) != set(assets):
        raise PreparationError("native Amber validation did not cover every advertised asset")
    return {
        "status": "pass", "canonical": canonical,
        "assets": validation_assets, "representative_numerical": numerical,
    }


def read_tagged_genesis_restart(path: Path) -> tuple[object, tuple[float, float, float]]:
    """Read coordinates and an orthorhombic box from a tagged GENESIS restart."""
    import numpy as np

    data = path.read_bytes()
    if data[:8] != b"GD150608":
        raise PreparationError("%s is not a tagged GENESIS restart" % path)
    records = {}
    offset = 8
    while offset < len(data):
        if offset + 8 > len(data):
            raise PreparationError("truncated GENESIS restart record in %s" % path)
        size = struct.unpack_from("=q", data, offset)[0]
        if size < 57 or offset + size > len(data):
            raise PreparationError("invalid GENESIS restart record in %s" % path)
        kind = chr(data[offset + 8])
        try:
            tag = data[offset + 9:offset + 49].decode("ascii").rstrip(" \0")
        except UnicodeDecodeError as error:
            raise PreparationError("invalid GENESIS restart tag in %s: %s" % (path, error))
        if not tag or tag in records:
            raise PreparationError("invalid or duplicate GENESIS restart tag %r" % tag)
        records[tag] = (kind, data[offset + 49:offset + size - 8])
        offset += size
    if offset != len(data):
        raise PreparationError("GENESIS restart record lengths do not cover %s" % path)

    try:
        natom = struct.unpack("=i", records["num_atoms"][1])[0]
    except (KeyError, struct.error) as error:
        raise PreparationError("GENESIS restart has no valid num_atoms record: %s" % error)

    def doubles(tag: str, count: int | None = None) -> tuple[float, ...]:
        try:
            kind, payload = records[tag]
        except KeyError:
            raise PreparationError("GENESIS restart is missing %s" % tag) from None
        if kind != "w" or len(payload) % 8:
            raise PreparationError("GENESIS restart %s is not a float64 record" % tag)
        values = struct.unpack("=%dd" % (len(payload) // 8), payload)
        if count is not None and len(values) != count:
            raise PreparationError("GENESIS restart %s has the wrong length" % tag)
        return values

    coordinates = np.asarray([
        doubles("coord_x", natom),
        doubles("coord_y", natom),
        doubles("coord_z", natom),
    ], dtype=float).T
    box = tuple(doubles(tag, 1)[0]
                for tag in ("box_size_x", "box_size_y", "box_size_z"))
    if any(not math.isfinite(value) or value <= 0.0 for value in box):
        raise PreparationError("GENESIS restart has an invalid periodic box")
    return coordinates, box


def is_charmm_water(residue: object) -> bool:
    return residue.name.upper() in WATER_NAMES and sorted(
        atom.atomic_number for atom in residue.atoms
    ) == [1, 1, 8]


def is_charmm_monatomic_ion(residue: object) -> bool:
    return len(residue.atoms) == 1 and residue.atoms[0].atomic_number in (11, 12, 17, 19, 20)


def select_native_apoa1_waters(
    source: object, count: int,
) -> tuple[list[tuple], dict]:
    """Select source-residue-stable ApoA1 TIP3 sites for native CHARMM SOD."""
    import numpy as np

    coordinates = np.asarray(source.coordinates, dtype=float)
    box = tuple(float(value) for value in source.box[:3])
    seed = "%s:apoa1" % NEUTRALIZATION_SEED_PREFIX
    waters = []
    existing_ions = []
    existing_ion_counts = {}
    solute_heavy = []
    candidates = []
    for residue in source.residues:
        if is_charmm_water(residue):
            waters.append(residue)
            oxygen = next(atom for atom in residue.atoms if atom.atomic_number == 8)
            rank = hashlib.sha256(
                seed.encode("ascii") + b"\0" + str(residue.idx).encode("ascii")
            ).digest()
            candidates.append(
                (rank, residue.idx, residue, oxygen, coordinates[oxygen.idx].copy())
            )
        elif is_charmm_monatomic_ion(residue):
            ion = residue.atoms[0]
            existing_ions.append(coordinates[ion.idx])
            existing_ion_counts[residue.name] = existing_ion_counts.get(residue.name, 0) + 1
        else:
            solute_heavy.extend(
                coordinates[atom.idx] for atom in residue.atoms if atom.atomic_number != 1
            )
    candidates.sort(key=lambda item: (item[0], item[1]))
    solute_heavy = np.asarray(solute_heavy, dtype=float)
    box_array = np.asarray(box, dtype=float)
    ion_cutoff2 = (MIN_ION_SEPARATION_NM * 10.0) ** 2
    solute_cutoff2 = (MIN_ION_SOLUTE_DISTANCE_NM * 10.0) ** 2
    selected = []
    for candidate in candidates:
        position = candidate[4]
        if any(pbc_distance_squared(position, other, box) < ion_cutoff2
               for other in existing_ions):
            continue
        if any(pbc_distance_squared(position, other[4], box) < ion_cutoff2
               for other in selected):
            continue
        delta = position - solute_heavy
        delta -= box_array * np.floor(delta / box_array + 0.5)
        if np.any(np.einsum("ij,ij->i", delta, delta) < solute_cutoff2):
            continue
        selected.append(candidate)
        if len(selected) == count:
            break
    if len(selected) != count:
        raise PreparationError("found only %d of %d safe ApoA1 ion sites" % (
            len(selected), count,
        ))
    metadata = {
        "method": "replace_waters",
        "seed": seed,
        "tolerance_e": FORMAL_CHARGE_TOLERANCE_E,
        "min_ion_separation_nm": MIN_ION_SEPARATION_NM,
        "min_ion_solute_distance_nm": MIN_ION_SOLUTE_DISTANCE_NM,
        "neutralization_ion": "Na+",
        "neutralization_ion_charge_e": 1,
        "neutralization_ion_count": count,
        "replaced_water_count": count,
        "existing_ions": dict(sorted(existing_ion_counts.items())),
        "water_count_before": len(waters),
        "water_count_after": len(waters) - count,
        "selected_water_residue_indices": [item[1] for item in selected],
        "selected_water_oxygen_positions_nm": [
            [float(value) / 10.0 for value in item[4]] for item in selected
        ],
    }
    return selected, metadata


def remove_charmm_donor_acceptor_terms(structure: object, atom_indices: set[int]) -> None:
    # Structure.strip removes bonded terms that contain deleted hydrogens.  PSF
    # donor/acceptor records are not bonded terms, so remove the selected water
    # records explicitly while leaving every unrelated diagnostic record intact.
    for name in ("donors", "acceptors"):
        items = getattr(structure, name)
        fields = CHARMM_TERM_FIELDS[name]
        for index in range(len(items) - 1, -1, -1):
            if any(getattr(items[index], field).idx in atom_indices for field in fields):
                del items[index]


def apply_native_apoa1_replacements(
    structure: object, selected_residue_indices: set[int],
) -> tuple[set[int], set[int], set[int]]:
    selected_atoms = {
        atom.idx for residue_index in selected_residue_indices
        for atom in structure.residues[residue_index].atoms
    }
    selected_oxygens = {
        next(atom.idx for atom in structure.residues[residue_index].atoms
             if atom.atomic_number == 8)
        for residue_index in selected_residue_indices
    }
    remove_charmm_donor_acceptor_terms(structure, selected_atoms)
    deleted_hydrogens = set()
    for residue_index in sorted(selected_residue_indices):
        residue = structure.residues[residue_index]
        if not is_charmm_water(residue):
            raise PreparationError("selected ApoA1 residue %d is not native TIP3" % residue_index)
        oxygen = next(atom for atom in residue.atoms if atom.atomic_number == 8)
        deleted_hydrogens.update(
            atom.idx for atom in residue.atoms if atom.atomic_number == 1
        )
        residue.name = "SOD"
        oxygen.name = "SOD"
        oxygen.type = 190
        oxygen.charge = 1.0
        oxygen.mass = 22.989770
        oxygen.atomic_number = 11
    structure.strip([atom.idx in deleted_hydrogens for atom in structure.atoms])
    return selected_atoms, selected_oxygens, deleted_hydrogens


def repartition_native_charmm_hmr(structure: object) -> dict:
    before = math.fsum(float(atom.mass) for atom in structure.atoms)
    changed_hydrogens = set()
    changed_heavy_atoms = set()
    for atom in structure.atoms:
        if atom.atomic_number != 1 or is_charmm_water(atom.residue):
            continue
        heavy = []
        for bond in atom.bonds:
            partner = bond.atom2 if bond.atom1 is atom else bond.atom1
            if partner.atomic_number != 1:
                heavy.append(partner)
        if len(heavy) != 1:
            raise PreparationError("ApoA1 hydrogen %d has %d bonded heavy atoms" % (
                atom.idx, len(heavy),
            ))
        delta = HMR_HYDROGEN_MASS_DA - float(atom.mass)
        if delta < -1.0e-9 or float(heavy[0].mass) <= delta:
            raise PreparationError("invalid ApoA1 HMR transfer at atom %d" % atom.idx)
        atom.mass = HMR_HYDROGEN_MASS_DA
        heavy[0].mass -= delta
        changed_hydrogens.add(atom.idx)
        changed_heavy_atoms.add(heavy[0].idx)
    after = math.fsum(float(atom.mass) for atom in structure.atoms)
    if not changed_hydrogens or abs(after - before) > 1.0e-6:
        raise PreparationError("native ApoA1 HMR failed mass conservation")
    changed = changed_hydrogens | changed_heavy_atoms
    return {
        "method": "native-topology heavy-atom mass transfer; rigid waters excluded",
        "hydrogen_mass_dalton": HMR_HYDROGEN_MASS_DA,
        "particles_with_changed_mass": len(changed),
        "water_particles_with_changed_mass": 0,
        "total_mass_delta_dalton": after - before,
    }


def mutate_native_apoa1(
    source_topology: Path, source_restart: Path,
) -> tuple[object, object, object, dict, dict, dict]:
    """Neutralize the pinned native CHARMM27 ApoA1 PSF without rebuilding it."""
    import parmed as pmd

    coordinates, box = read_tagged_genesis_restart(source_restart)
    source = pmd.load_file(str(source_topology))
    if len(source.atoms) != len(coordinates):
        raise PreparationError("native ApoA1 PSF/restart atom count mismatch")
    source.coordinates = coordinates
    source.box = [*box, 90.0, 90.0, 90.0]
    source_charge = math.fsum(float(atom.charge) for atom in source.atoms)
    source_formal = integer_charge(source_charge, "native ApoA1 source")
    if source_formal != -14:
        raise PreparationError("native ApoA1 formal charge changed: expected -14, got %d" % (
            source_formal,
        ))
    selected, neutralization = select_native_apoa1_waters(source, -source_formal)
    selected_residue_indices = {item[1] for item in selected}

    # Load twice from disk.  ParmEd deepcopy does not preserve the source PSF's
    # per-atom CHEQ property strings, which are part of the native topology.
    neutral = pmd.load_file(str(source_topology))
    neutral.coordinates = coordinates
    neutral.box = [*box, 90.0, 90.0, 90.0]
    selected_atoms, selected_oxygens, deleted_hydrogens = apply_native_apoa1_replacements(
        neutral, selected_residue_indices,
    )
    hmr = pmd.load_file(str(source_topology))
    hmr.coordinates = coordinates
    hmr.box = [*box, 90.0, 90.0, 90.0]
    apply_native_apoa1_replacements(hmr, selected_residue_indices)
    hmr_metadata = repartition_native_charmm_hmr(hmr)

    post_charge = math.fsum(float(atom.charge) for atom in neutral.atoms)
    if integer_charge(post_charge, "neutral native ApoA1") != 0:
        raise PreparationError("neutralized native ApoA1 has nonzero formal charge")
    if abs(post_charge) > POST_REPLACEMENT_TOLERANCE_E:
        raise PreparationError("neutralized native ApoA1 charge is %.12g e" % post_charge)
    neutralization.update({
        "source_charge_e": source_charge,
        "source_formal_charge_e": source_formal,
        "pre_neutralization_charge_e": source_charge,
        "pre_neutralization_formal_charge_e": source_formal,
        "post_neutralization_charge_e": post_charge,
        "post_neutralization_formal_charge_e": 0,
    })
    preservation = {
        "source": source,
        "source_atom_count": len(source.atoms),
        "selected_source_atoms": selected_atoms,
        "selected_source_oxygens": selected_oxygens,
        "deleted_hydrogen_count": len(deleted_hydrogens),
    }
    return source, neutral, hmr, neutralization, hmr_metadata, preservation


def write_namd_binary_coordinates(path: Path, coordinates: object) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack("=i", len(coordinates)))
        for xyz in coordinates:
            stream.write(struct.pack("=3d", *map(float, xyz)))


def read_namd_binary_coordinates(path: Path) -> object:
    import numpy as np

    data = path.read_bytes()
    if len(data) < 4:
        raise PreparationError("truncated NAMD binary coordinate file: %s" % path)
    natom = struct.unpack_from("=i", data)[0]
    if natom < 1 or len(data) != 4 + natom * 24:
        raise PreparationError("invalid NAMD binary coordinate length: %s" % path)
    return np.asarray(struct.unpack_from("=%dd" % (natom * 3), data, 4)).reshape(natom, 3)


def read_pdb_scaffold(path: Path) -> tuple[list[tuple[str, str, str]], list[list[float]]]:
    identities = []
    coordinates = []
    for line in path.read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise PreparationError("short atom record in %s" % path)
        try:
            coordinates.append([
                float(line[30:38]), float(line[38:46]), float(line[46:54]),
            ])
        except ValueError as error:
            raise PreparationError("invalid coordinate in %s: %s" % (path, error))
        identities.append((line[12:16].strip(), line[17:21].strip(), line[72:76].strip()))
    if not identities:
        raise PreparationError("PDB scaffold has no atoms: %s" % path)
    return identities, coordinates


def parameterize_charmm_psf(structure: object, parameters: object) -> object:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure.load_parameters(parameters)
    structure.flags = ["CMAP", "XPLOR"]
    return structure


def export_native_apoa1(
    directory: Path, neutral: object, hmr: object,
    source_rtf: Path, source_prm: Path,
) -> tuple[dict, object]:
    import parmed as pmd
    from parmed.charmm import CharmmParameterSet, CharmmPsfFile

    rtf_name = source_rtf.name
    prm_name = source_prm.name
    shutil.copy2(source_rtf, directory / rtf_name)
    shutil.copy2(source_prm, directory / prm_name)
    # The legacy RTF contains a DELETE ACCE directive used only while building
    # patches.  We preserve the already patched PSF, so ParmEd's informational
    # warning that its parameter reader does not apply this build directive is
    # irrelevant to direct-asset parameter coverage.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parameters = CharmmParameterSet(str(directory / rtf_name), str(directory / prm_name))

    neutral.flags = ["CMAP"]
    hmr.flags = ["CMAP"]
    CharmmPsfFile.write_psf(neutral, str(directory / "system.psf"))
    temporary_hmr = directory / ".system_hmr_numeric.psf"
    CharmmPsfFile.write_psf(hmr, str(temporary_hmr))

    numeric = pmd.load_file(str(directory / "system.psf"))
    numeric.coordinates = neutral.coordinates
    numeric.box = list(neutral.box)
    numeric_hmr = pmd.load_file(str(temporary_hmr))
    numeric_hmr.coordinates = neutral.coordinates
    numeric_hmr.box = list(neutral.box)
    parameterize_charmm_psf(numeric, parameters)
    parameterize_charmm_psf(numeric_hmr, parameters)
    CharmmPsfFile.write_psf(numeric, str(directory / "system_xplor.psf"))
    CharmmPsfFile.write_psf(numeric_hmr, str(directory / "system_hmr_xplor.psf"))
    temporary_hmr.unlink()

    pmd.charmm.CharmmCrdFile.write(neutral, str(directory / "system.crd"))
    write_namd_binary_coordinates(directory / "system.coor", neutral.coordinates)
    # NAMD requires a text `coordinates` file even when authoritative binary
    # coordinates follow.  Atom order and CHARMM SEGIDs are retained; the PDB's
    # three-decimal positions are only a scaffold and are never authoritative.
    pmd.formats.PDBFile.write(
        neutral, str(directory / "system.pdb"), renumber=True, charmm=True,
        use_hetatoms=False, increase_tercount=False, write_links=False,
    )

    assets = {
        "GENESIS": {
            "format": "CHARMM",
            "topology": "system.psf",
            "coordinates": "system.crd",
            "topology_definitions": [rtf_name],
            "parameters": [prm_name],
        },
        "NAMD": {
            "format": "CHARMM",
            "topology": "system_xplor.psf",
            "topology_hmr": "system_hmr_xplor.psf",
            "coordinates": "system.coor",
            "coordinate_reference": "system.pdb",
            "parameters": [prm_name],
        },
    }
    return assets, parameters


def charmm_term_indices(items: object, fields: tuple[str, ...]) -> list[tuple[int, ...]]:
    return [tuple(getattr(item, field).idx for field in fields) for item in items]


def charmm_atom_identity(atom: object) -> tuple:
    return (
        atom.name, str(atom.type), atom.residue.name, atom.residue.number,
        atom.residue.segid, float(atom.charge), tuple(getattr(atom, "props", ())),
    )


def validate_apoa1_source_preservation(
    source: object, output: object, preservation: dict,
) -> dict:
    selected_atoms = preservation["selected_source_atoms"]
    selected_oxygens = preservation["selected_source_oxygens"]
    expected = [
        atom for atom in source.atoms
        if atom.idx not in selected_atoms or atom.idx in selected_oxygens
    ]
    if len(expected) != len(output.atoms):
        raise PreparationError("ApoA1 output atom count does not match selected waters")
    output_to_source = {}
    converted = 0
    for observed, original in zip(output.atoms, expected):
        output_to_source[observed.idx] = original.idx
        if original.idx in selected_oxygens:
            if not (
                observed.name == "SOD"
                and str(observed.type) == "190"
                and observed.residue.name == "SOD"
                and observed.residue.number == original.residue.number
                and observed.residue.segid == original.residue.segid
                and abs(float(observed.charge) - 1.0) < 1.0e-12
                and abs(float(observed.mass) - 22.989770) <= 5.1e-5
                and observed.atomic_number == 11
                and tuple(getattr(observed, "props", ()))
                    == tuple(getattr(original, "props", ()))
            ):
                raise PreparationError("selected ApoA1 water oxygen is not canonical SOD")
            converted += 1
            continue
        expected_identity = (
            original.name, str(original.type), original.residue.name,
            original.residue.number, original.residue.segid, float(original.charge),
            tuple(getattr(original, "props", ())),
        )
        if (charmm_atom_identity(observed) != expected_identity
                or abs(float(observed.mass) - float(original.mass)) > 1.0e-12
                or observed.atomic_number != original.atomic_number):
            raise PreparationError(
                "native ApoA1 atom parameters changed outside selected waters at %d" % original.idx
            )
    if converted != len(selected_oxygens):
        raise PreparationError("ApoA1 neutralization converted the wrong number of waters")

    term_report = {}
    for name, fields in CHARMM_TERM_FIELDS.items():
        original_terms = Counter(
            tuple(getattr(term, field).idx for field in fields)
            for term in getattr(source, name)
            if not any(getattr(term, field).idx in selected_atoms for field in fields)
        )
        output_terms = Counter(
            tuple(output_to_source[getattr(term, field).idx] for field in fields)
            for term in getattr(output, name)
        )
        if output_terms != original_terms:
            raise PreparationError("native ApoA1 %s terms changed outside selected waters" % name)
        term_report[name] = {
            "source": len(getattr(source, name)),
            "retained": len(getattr(output, name)),
        }
    removed_atoms = selected_atoms - selected_oxygens
    original_groups = Counter(
        (group.atom.idx, group.type, group.move)
        for group in source.groups if group.atom.idx not in removed_atoms
    )
    output_groups = Counter(
        (output_to_source[group.atom.idx], group.type, group.move) for group in output.groups
    )
    if output_groups != original_groups:
        raise PreparationError("native ApoA1 charge groups changed")
    return term_report


def validate_apoa1_xplor_mapping(
    numeric: object, xplor: object, parameters: object,
) -> None:
    if len(numeric.atoms) != len(xplor.atoms):
        raise PreparationError("numeric/XPLOR ApoA1 PSF atom counts differ")
    for normal, named in zip(numeric.atoms, xplor.atoms):
        try:
            expected_type = parameters.atom_types_int[int(normal.type)].name
        except (KeyError, TypeError, ValueError) as error:
            raise PreparationError("cannot map ApoA1 numeric CHARMM type: %s" % error)
        if (
            normal.name != named.name
            or normal.residue.name != named.residue.name
            or normal.residue.number != named.residue.number
            or normal.residue.segid != named.residue.segid
            or expected_type != str(named.type)
            or abs(float(normal.charge) - float(named.charge)) > 1.0e-12
            or abs(float(normal.mass) - float(named.mass)) > 1.0e-12
            or normal.atomic_number != named.atomic_number
            or tuple(getattr(normal, "props", ())) != tuple(getattr(named, "props", ()))
        ):
            raise PreparationError("numeric/XPLOR ApoA1 atom metadata differs")
    for name, fields in CHARMM_TERM_FIELDS.items():
        if charmm_term_indices(getattr(numeric, name), fields) != charmm_term_indices(
                getattr(xplor, name), fields):
            raise PreparationError("numeric/XPLOR ApoA1 %s terms differ" % name)
    if ([(group.atom.idx, group.type, group.move) for group in numeric.groups]
            != [(group.atom.idx, group.type, group.move) for group in xplor.groups]):
        raise PreparationError("numeric/XPLOR ApoA1 charge groups differ")


def validate_apoa1_hmr(normal: object, hmr: object, hmr_metadata: dict) -> tuple[float, float]:
    if len(normal.atoms) != len(hmr.atoms):
        raise PreparationError("ApoA1 HMR PSF changed atom count")
    changed = 0
    water_changes = 0
    for ordinary, repartitioned in zip(normal.atoms, hmr.atoms):
        # A PSF has no explicit element column.  Do not compare ParmEd's element
        # inference after HMR: a correctly reduced heavy-atom mass can resemble
        # another element.  All serialized non-mass identity fields are checked.
        if charmm_atom_identity(ordinary) != charmm_atom_identity(repartitioned):
            raise PreparationError("ApoA1 HMR changed non-mass atom metadata")
        expected_mass = float(ordinary.mass)
        if ordinary.atomic_number == 1 and not is_charmm_water(ordinary.residue):
            expected_mass = HMR_HYDROGEN_MASS_DA
        elif ordinary.atomic_number != 1:
            expected_mass -= math.fsum(
                HMR_HYDROGEN_MASS_DA - float(partner.mass)
                for bond in ordinary.bonds
                for partner in (bond.atom2 if bond.atom1 is ordinary else bond.atom1,)
                if partner.atomic_number == 1 and not is_charmm_water(partner.residue)
            )
        if abs(float(repartitioned.mass) - expected_mass) > 5.1e-5:
            raise PreparationError("unexpected ApoA1 HMR mass at atom %d" % ordinary.idx)
        if abs(float(repartitioned.mass) - float(ordinary.mass)) > 5.1e-5:
            changed += 1
            if is_charmm_water(ordinary.residue):
                water_changes += 1
    if changed != hmr_metadata["particles_with_changed_mass"] or water_changes != 0:
        raise PreparationError("ApoA1 HMR changed-particle metadata is inconsistent")
    for name, fields in CHARMM_TERM_FIELDS.items():
        if charmm_term_indices(getattr(normal, name), fields) != charmm_term_indices(
                getattr(hmr, name), fields):
            raise PreparationError("ApoA1 HMR changed %s terms" % name)
    if ([(group.atom.idx, group.type, group.move) for group in normal.groups]
            != [(group.atom.idx, group.type, group.move) for group in hmr.groups]):
        raise PreparationError("ApoA1 HMR changed charge groups")
    normal_mass = math.fsum(float(atom.mass) for atom in normal.atoms)
    hmr_mass = math.fsum(float(atom.mass) for atom in hmr.atoms)
    if abs(hmr_mass - normal_mass) > 1.0e-6:
        raise PreparationError("serialized ApoA1 HMR PSF changed total mass")
    if abs(hmr_mass - normal_mass - hmr_metadata["total_mass_delta_dalton"]) > 1.0e-6:
        raise PreparationError("ApoA1 HMR mass metadata does not match the direct asset")
    return normal_mass, hmr_mass


def charmm_exception_count(structure: object) -> int:
    pairs = {
        tuple(sorted((bond.atom1.idx, bond.atom2.idx))) for bond in structure.bonds
    }
    pairs.update(
        tuple(sorted((angle.atom1.idx, angle.atom3.idx))) for angle in structure.angles
    )
    pairs.update(
        tuple(sorted((dihedral.atom1.idx, dihedral.atom4.idx)))
        for dihedral in structure.dihedrals
    )
    return len(pairs)


def validate_native_apoa1(
    directory: Path,
    source: object,
    neutral: object,
    hmr_metadata: dict,
    preservation: dict,
    parameters: object,
) -> dict:
    import parmed as pmd

    numeric = pmd.load_file(str(directory / "system.psf"))
    xplor = pmd.load_file(str(directory / "system_xplor.psf"))
    xplor_hmr = pmd.load_file(str(directory / "system_hmr_xplor.psf"))
    crd = pmd.charmm.CharmmCrdFile(str(directory / "system.crd"))
    numeric.coordinates = crd.coordinates
    numeric.box = list(neutral.box)
    xplor.coordinates = crd.coordinates
    xplor.box = list(neutral.box)
    xplor_hmr.box = list(neutral.box)

    if set(numeric.flags) != {"CHEQ", "CMAP"}:
        raise PreparationError("numeric ApoA1 PSF lost CHEQ/CMAP flags")
    for label, structure in (("normal", xplor), ("HMR", xplor_hmr)):
        if set(structure.flags) != {"CHEQ", "CMAP", "XPLOR"}:
            raise PreparationError("%s ApoA1 XPLOR PSF lost native flags" % label)
    preservation_terms = validate_apoa1_source_preservation(source, numeric, preservation)
    validate_apoa1_xplor_mapping(numeric, xplor, parameters)
    normal_mass, hmr_mass = validate_apoa1_hmr(xplor, xplor_hmr, hmr_metadata)

    raw_charge = math.fsum(float(atom.charge) for atom in numeric.atoms)
    if integer_charge(raw_charge, "reloaded ApoA1 PSF") != 0:
        raise PreparationError("reloaded ApoA1 PSF is not neutral")
    reference_coordinates = neutral.coordinates
    try:
        crd_delta = coordinate_delta(reference_coordinates, numeric.coordinates)
    except PreparationError as error:
        raise PreparationError("ApoA1 EXT CRD validation failed: %s" % error)
    namd_coordinates = read_namd_binary_coordinates(directory / "system.coor")
    try:
        namd_delta = coordinate_delta(reference_coordinates, namd_coordinates)
    except PreparationError as error:
        raise PreparationError("ApoA1 NAMD binary validation failed: %s" % error)
    pdb_identities, pdb_coordinates = read_pdb_scaffold(directory / "system.pdb")
    expected_identities = [
        (atom.name, atom.residue.name, atom.residue.segid) for atom in numeric.atoms
    ]
    if pdb_identities != expected_identities:
        raise PreparationError("NAMD PDB scaffold changed ApoA1 atom order/identity")
    try:
        pdb_delta = coordinate_delta(reference_coordinates, pdb_coordinates)
    except PreparationError as error:
        raise PreparationError("ApoA1 PDB scaffold validation failed: %s" % error)
    if crd_delta > 5.1e-10 or namd_delta > 5.1e-10 or pdb_delta > 5.1e-4:
        raise PreparationError("ApoA1 coordinate export exceeded its file precision")
    box = [float(value) for value in neutral.box[:3]]
    if any(abs(float(numeric.box[index]) - box[index]) > 1.0e-10 for index in range(3)):
        raise PreparationError("ApoA1 periodic box changed during reload")
    if len(numeric.atoms) != len(neutral.atoms):
        raise PreparationError("ApoA1 direct asset changed atom count")

    canonical = {
        "atom_count": len(numeric.atoms),
        "raw_charge_e": raw_charge,
        "formal_charge_e": integer_charge(raw_charge, "native ApoA1 canonical asset"),
        "box_angstrom": box,
        "bonded_term_counts": {
            "bonds": len(numeric.bonds),
            "angles": len(numeric.angles),
            "proper_dihedrals": len(numeric.dihedrals),
            "impropers": len(numeric.impropers),
            "cmaps": len(numeric.cmaps),
        },
        "exception_count": charmm_exception_count(numeric),
        "constraint_count": 0,
    }
    genesis = asset_record(
        len(numeric.atoms), raw_charge, box, crd_delta, normal_mass, None,
    )
    namd = asset_record(
        len(numeric.atoms), raw_charge, box, max(namd_delta, pdb_delta),
        normal_mass, hmr_mass,
    )
    # Numeric/XPLOR atom-type mapping and static HMR are direct topology
    # validations, not a force-field conversion.  Cross-engine numerical
    # equivalence is therefore neither required nor claimed here.
    numerical = {
        "status": "not_run",
        "reference_variant": None,
        "energy_delta_kj_mol": None,
        "force_relative_rms": None,
        "cutoff_nm": None,
        "ewald_error_tolerance": None,
        "platform": None,
    }
    return {
        "status": "pass",
        "canonical": canonical,
        "assets": {"GENESIS": genesis, "NAMD": namd},
        "representative_numerical": numerical,
    }, preservation_terms


def parse_surgical_charmm_psf(path: Path) -> dict:
    """Parse the PSF fields needed for byte-preserving C36 transformations."""
    raw = path.read_bytes()
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as error:
        raise PreparationError("non-ASCII CHARMM PSF %s: %s" % (path, error))
    lines = text.splitlines(keepends=True)
    try:
        atom_header = next(index for index, line in enumerate(lines) if "!NATOM" in line)
        atom_count = int(lines[atom_header].split()[0])
    except (StopIteration, ValueError, IndexError) as error:
        raise PreparationError("invalid CHARMM PSF atom header in %s: %s" % (path, error))
    atom_start = atom_header + 1
    atom_end = atom_start + atom_count
    if atom_end > len(lines):
        raise PreparationError("truncated CHARMM PSF atom table in %s" % path)

    atoms = []
    residue_ordinal = 0
    previous_residue = None
    for expected, line in enumerate(lines[atom_start:atom_end], 1):
        tokens = list(re.finditer(r"\S+", line))
        fields = [token.group() for token in tokens]
        if len(fields) < 8:
            raise PreparationError("malformed CHARMM PSF atom %d in %s" % (expected, path))
        try:
            index = int(fields[0])
            charge = Decimal(fields[6])
            mass = Decimal(fields[7])
        except (ValueError, ArithmeticError) as error:
            raise PreparationError("invalid CHARMM PSF atom %d in %s: %s" % (
                expected, path, error,
            ))
        if index != expected:
            raise PreparationError("nonsequential CHARMM PSF atom %d in %s" % (expected, path))
        residue_key = (fields[1], fields[2])
        if residue_key != previous_residue:
            residue_ordinal += 1
            previous_residue = residue_key
        atoms.append({
            "index": index,
            "segid": fields[1],
            "resid": fields[2],
            "resname": fields[3],
            "atomname": fields[4],
            "atomtype": fields[5],
            "charge": charge,
            "mass": mass,
            "mass_span": tokens[7].span(),
            "next_field_start": (
                tokens[8].span()[0] if len(tokens) > 8 else len(line.rstrip("\r\n"))
            ),
            "line": line,
            "residue_ordinal": residue_ordinal,
        })
    return {
        "raw": raw,
        "lines": lines,
        "flags": tuple(lines[0].split()[1:]),
        "atom_header": atom_header,
        "atom_start": atom_start,
        "atom_end": atom_end,
        "atoms": atoms,
    }


def charmm_psf_integer_terms(parsed: dict, marker: str, width: int) -> list[tuple[int, ...]]:
    lines = parsed["lines"]
    try:
        header = next(index for index, line in enumerate(lines) if marker in line)
        count = int(lines[header].split()[0])
    except (StopIteration, ValueError, IndexError) as error:
        raise PreparationError("invalid %s section in CHARMM PSF: %s" % (marker, error))
    values = []
    line_index = header + 1
    expected = count * width
    while len(values) < expected and line_index < len(lines):
        try:
            values.extend(int(value) for value in lines[line_index].split())
        except ValueError as error:
            raise PreparationError("invalid integer in CHARMM PSF %s section: %s" % (
                marker, error,
            ))
        line_index += 1
    if len(values) != expected:
        raise PreparationError("CHARMM PSF %s section has inconsistent length" % marker)
    return [tuple(values[index:index + width]) for index in range(0, expected, width)]


def charmm_psf_section_count(parsed: dict, marker: str) -> int:
    try:
        return int(next(line for line in parsed["lines"] if marker in line).split()[0])
    except (StopIteration, ValueError, IndexError) as error:
        raise PreparationError("invalid %s count in CHARMM PSF: %s" % (marker, error))


def native_charmm_psf_summary(parsed: dict) -> dict:
    atoms = parsed["atoms"]
    bonds = charmm_psf_integer_terms(parsed, "!NBOND", 2)
    angles = charmm_psf_integer_terms(parsed, "!NTHETA", 3)
    dihedrals = charmm_psf_integer_terms(parsed, "!NPHI", 4)
    exceptions = {tuple(sorted(term)) for term in bonds}
    exceptions.update(tuple(sorted((term[0], term[2]))) for term in angles)
    exceptions.update(tuple(sorted((term[0], term[3]))) for term in dihedrals)

    residues = []
    previous = None
    for atom in atoms:
        key = (atom["segid"], atom["resid"])
        if key != previous:
            residues.append([key, atom["resname"], 0])
            previous = key
        elif residues[-1][1] != atom["resname"]:
            raise PreparationError("CHARMM PSF residue name changes within one residue")
        residues[-1][2] += 1
    water_count = sum(
        count == 3 and name.upper() in WATER_NAMES for _, name, count in residues
    )
    ion_names = frozenset(("LIT", "SOD", "MG", "POT", "CAL", "RUB", "CES", "BAR", "ZN2", "CD2", "CLA"))
    existing_ions = Counter(
        name for _, name, count in residues if count == 1 and name.upper() in ion_names
    )
    return {
        "atom_count": len(atoms),
        "raw_charge": sum((atom["charge"] for atom in atoms), Decimal(0)),
        "normal_mass": sum((atom["mass"] for atom in atoms), Decimal(0)),
        "water_count": water_count,
        "existing_ions": dict(sorted(existing_ions.items())),
        "bonded_term_counts": {
            "bonds": len(bonds),
            "angles": len(angles),
            "proper_dihedrals": len(dihedrals),
            "impropers": charmm_psf_section_count(parsed, "!NIMPHI"),
            "cmaps": charmm_psf_section_count(parsed, "!NCRTERM"),
        },
        "exception_count": len(exceptions),
        "pdb_identities": [
            (atom["atomname"], atom["resname"], atom["segid"]) for atom in atoms
        ],
    }


def charmm_mass_type_map(paths: tuple[Path, ...]) -> dict[int, str]:
    mapping = {}
    for path in paths:
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            match = re.match(r"^\s*MASS\s+(-?\d+)\s+(\S+)\s+([0-9.Ee+-]+)", line, re.I)
            if not match:
                continue
            number = int(match.group(1))
            name = match.group(2)
            if number in mapping and mapping[number] != name:
                raise PreparationError(
                    "numeric CHARMM type %d maps to both %s and %s" % (
                        number, mapping[number], name,
                    )
                )
            mapping[number] = name
    return mapping


def charmm_named_types(paths: tuple[Path, ...]) -> set[str]:
    names = set()
    for path in paths:
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            match = re.match(r"^\s*MASS\s+-?\d+\s+(\S+)\s+([0-9.Ee+-]+)", line, re.I)
            if match:
                names.add(match.group(1))
    return names


def write_native_charmm_xplor_psf(
    source: Path, destination: Path, topology_cards: tuple[Path, ...],
) -> None:
    parsed = parse_surgical_charmm_psf(source)
    if "XPLOR" in parsed["flags"]:
        shutil.copyfile(source, destination)
        return
    if "EXT" in parsed["flags"]:
        raise PreparationError("numeric EXT PSF conversion is not implemented: %s" % source)
    mapping = charmm_mass_type_map(topology_cards)
    converted = list(parsed["lines"])
    converted[0] = parsed["lines"][0].rstrip("\r\n") + " XPLOR\n"
    for atom, line_index in zip(
            parsed["atoms"], range(parsed["atom_start"], parsed["atom_end"])):
        try:
            expected = mapping[int(atom["atomtype"])]
        except (KeyError, ValueError) as error:
            raise PreparationError("cannot map native numeric CHARMM type %s: %s" % (
                atom["atomtype"], error,
            ))
        line = parsed["lines"][line_index]
        if len(expected) > 6 or line[29:35].strip() != atom["atomtype"]:
            raise PreparationError("unsupported standard CHARMM PSF atom layout at %d" % (
                atom["index"],
            ))
        converted[line_index] = line[:29] + ("%-6s" % expected) + line[35:]
    destination.write_text("".join(converted), encoding="ascii")


def native_charmm_is_hydrogen(atom: dict) -> bool:
    return (
        Decimal("0.5") < atom["mass"] < Decimal("2.0")
        and (atom["atomname"].upper().startswith("H")
             or atom["atomtype"].upper().startswith("H"))
    )


def write_native_charmm_hmr_psf(source: Path, destination: Path) -> dict:
    """Surgically repartition only non-water hydrogens in a named-type PSF."""
    parsed = parse_surgical_charmm_psf(source)
    if "XPLOR" not in parsed["flags"]:
        raise PreparationError("NAMD HMR requires a named-type XPLOR PSF")
    atoms = parsed["atoms"]
    bonds = charmm_psf_integer_terms(parsed, "!NBOND", 2)
    adjacency = defaultdict(list)
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)

    water_indices = {
        atom["index"] for atom in atoms if atom["resname"].upper() == "TIP3"
    }
    hydrogens = {atom["index"] for atom in atoms if native_charmm_is_hydrogen(atom)}
    solute_hydrogens = sorted(hydrogens - water_indices)
    donor_hydrogens = defaultdict(list)
    for hydrogen_index in solute_hydrogens:
        heavy = [
            other for other in adjacency[hydrogen_index]
            if other not in hydrogens and other not in water_indices
        ]
        if len(heavy) != 1:
            atom = atoms[hydrogen_index - 1]
            raise PreparationError(
                "solute hydrogen %d %s:%s:%s has %d non-water donors" % (
                    hydrogen_index, atom["segid"], atom["resid"], atom["atomname"],
                    len(heavy),
                )
            )
        donor_hydrogens[heavy[0]].append(hydrogen_index)

    target = Decimal("3.02400")
    new_masses = [atom["mass"] for atom in atoms]
    for donor_index in sorted(donor_hydrogens):
        children = donor_hydrogens[donor_index]
        delta = sum(
            (target - atoms[index - 1]["mass"] for index in children), Decimal(0)
        )
        new_donor_mass = atoms[donor_index - 1]["mass"] - delta
        if new_donor_mass <= 0:
            raise PreparationError("HMR makes donor %d non-positive" % donor_index)
        new_masses[donor_index - 1] = new_donor_mass
        for index in children:
            new_masses[index - 1] = target

    changed_indices = set(solute_hydrogens) | set(donor_hydrogens)
    output_lines = list(parsed["lines"])
    for atom, new_mass in zip(atoms, new_masses):
        if atom["index"] not in changed_indices:
            continue
        value = format(new_mass, ".5f")
        start, end = atom["mass_span"]
        field_end = atom["next_field_start"]
        width = field_end - start
        if len(value) >= width:
            raise PreparationError("HMR mass does not fit atom %d field" % atom["index"])
        line = atom["line"]
        output_lines[parsed["atom_start"] + atom["index"] - 1] = (
            line[:start] + value + " " * (field_end - start - len(value)) + line[field_end:]
        )
    destination.write_text("".join(output_lines), encoding="ascii")

    output = parse_surgical_charmm_psf(destination)
    if len(output["atoms"]) != len(atoms):
        raise PreparationError("HMR changed CHARMM PSF atom count")
    if parsed["lines"][:parsed["atom_start"]] != output["lines"][:output["atom_start"]]:
        raise PreparationError("HMR changed CHARMM PSF header")
    if parsed["lines"][parsed["atom_end"]:] != output["lines"][output["atom_end"]:]:
        raise PreparationError("HMR changed CHARMM PSF topology sections")
    water_changes = 0
    for before, after, expected_mass in zip(atoms, output["atoms"], new_masses):
        for field in (
                "index", "segid", "resid", "resname", "atomname", "atomtype", "charge"):
            if before[field] != after[field]:
                raise PreparationError("HMR changed atom %d field %s" % (
                    before["index"], field,
                ))
        if after["mass"] != expected_mass:
            raise PreparationError("HMR mass did not serialize exactly at atom %d" % before["index"])
        if before["index"] in water_indices and after["mass"] != before["mass"]:
            water_changes += 1
        if before["index"] not in changed_indices and after["mass"] != before["mass"]:
            raise PreparationError("HMR changed unselected atom %d" % before["index"])
    source_mass = sum((atom["mass"] for atom in atoms), Decimal(0))
    output_mass = sum((atom["mass"] for atom in output["atoms"]), Decimal(0))
    if source_mass != output_mass or water_changes:
        raise PreparationError("HMR changed total mass or rigid-water masses")
    if sum((atom["charge"] for atom in atoms), Decimal(0)) != sum(
            (atom["charge"] for atom in output["atoms"]), Decimal(0)):
        raise PreparationError("HMR changed total charge")
    return {
        "hydrogen_mass_dalton": HMR_HYDROGEN_MASS_DA,
        "method": "native-topology heavy-atom mass transfer; rigid waters excluded",
        "particles_with_changed_mass": len(changed_indices),
        "water_particles_with_changed_mass": water_changes,
        "total_mass_delta_dalton": float(output_mass - source_mass),
        "solute_hydrogen_count": len(solute_hydrogens),
        "donor_count": len(donor_hydrogens),
        "normal_total_mass": source_mass,
        "hmr_total_mass": output_mass,
    }


def embedded_charmm_sections(path: Path) -> tuple[str, str, str]:
    lines = path.read_text(encoding="utf-8-sig").splitlines(keepends=True)
    sections = []
    for index, line in enumerate(lines):
        words = line.lower().split()
        if (len(words) < 3 or words[0] != "read" or words[1] not in ("rtf", "para")
                or words[2] != "card"):
            continue
        end = index + 1
        while end < len(lines) and lines[end].strip().upper() != "END":
            end += 1
        if end == len(lines):
            raise PreparationError("unterminated embedded %s card in %s" % (words[1], path))
        sections.append((words[1], "".join(lines[index + 1:end + 1])))
    if [kind for kind, _ in sections] != ["rtf", "para", "para"]:
        raise PreparationError("%s must contain RTF, parameter, and NBFIX cards" % path)
    return tuple(content for _, content in sections)


def charmm_nbfix_entries(content: str) -> list[tuple[str, str, str]]:
    entries = []
    active = False
    for line in content.splitlines(keepends=True):
        body = line.split("!", 1)[0].strip()
        if not body:
            continue
        if body.upper() == "NBFIX":
            active = True
            continue
        if body.upper() == "END":
            active = False
            continue
        if active:
            fields = body.split()
            if len(fields) < 4:
                raise PreparationError("malformed flattened CHARMM NBFIX entry: %s" % body)
            entries.append((fields[0], fields[1], line))
    return entries


def flatten_native_charmm_stream(
    system: str,
    stream: Path,
    topology_cards: tuple[Path, ...],
    directory: Path,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    topology_names = []
    parameter_names = []
    for source in topology_cards:
        destination = directory / source.name
        if destination.exists():
            raise PreparationError("duplicate flattened CHARMM topology name: %s" % source.name)
        shutil.copy2(source, destination)
        topology_names.append(source.name)

    rtf_text, parameter_text, nbfix_text = embedded_charmm_sections(stream)
    water_rtf = directory / "water_ions.rtf"
    water_parameter = directory / "water_ions.prm"
    water_nbfix = directory / "water_ions_nbfix.prm"
    water_rtf.write_text(rtf_text, encoding="utf-8")
    water_parameter.write_text(parameter_text, encoding="utf-8")
    topology_names.append(water_rtf.name)
    defined_types = charmm_named_types(tuple(directory / name for name in topology_names))

    entries = charmm_nbfix_entries(nbfix_text)
    retained = [line for first, second, line in entries
                if first in defined_types and second in defined_types]
    if not retained:
        raise PreparationError("flattened CHARMM stream has no applicable NBFIX entries")
    if len(retained) == len(entries):
        filtered_nbfix = nbfix_text
    else:
        filtered_nbfix = (
            "* NBFix terms from toppar_water_ions.str whose atom types are defined by the\n"
            "* DPPC native lipid plus water/ion parameter cards.\n"
            "*\n\nNBFIX\n" + "".join(retained) + "END\n"
        )
    water_nbfix.write_text(filtered_nbfix, encoding="utf-8")
    for first, second, _ in charmm_nbfix_entries(filtered_nbfix):
        if first not in defined_types or second not in defined_types:
            raise PreparationError("flattened NBFIX references an undefined atom type")
    parameter_names.extend((water_parameter.name, water_nbfix.name))
    return tuple(topology_names), tuple(parameter_names)


def write_native_charmm_crd(path: Path, atoms: list[dict], coordinates: object) -> None:
    if len(atoms) != len(coordinates):
        raise PreparationError("CHARMM coordinate atom counts differ")
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("* PSF-order coordinates exported from a GENESIS restart\n")
        stream.write("* Floating-point values use CHARMM EXT F20.10 fields\n")
        stream.write("*\n")
        stream.write("%10d  EXT\n" % len(atoms))
        for atom, xyz in zip(atoms, coordinates):
            integer_fields = (int(atom["index"]), int(atom["residue_ordinal"]))
            text_fields = tuple(str(atom[name]) for name in (
                "resname", "atomname", "segid",
            ))
            try:
                residue_id = int(atom["resid"])
            except (TypeError, ValueError) as error:
                raise PreparationError("CHARMM EXT residue ID must be an integer: %s" % error)
            if any(len(str(value)) > 10 for value in integer_fields):
                raise PreparationError("CHARMM EXT atom/residue ordinal exceeds I10")
            if any(not value or len(value) > 8 for value in text_fields):
                raise PreparationError("CHARMM EXT name is empty or exceeds A8")
            if len(str(residue_id)) > 8:
                raise PreparationError("CHARMM EXT residue ID exceeds I8")
            coordinate_fields = tuple("%20.10f" % float(value) for value in xyz)
            if (len(coordinate_fields) != 3
                    or any(len(value) != 20 or "nan" in value.lower()
                           or "inf" in value.lower() for value in coordinate_fields)):
                raise PreparationError("CHARMM EXT coordinate is non-finite or exceeds F20.10")
            line = (
                "%10d%10d  %-8s  %-8s%s%s%s  %-8s  %8d" % (
                    *integer_fields, text_fields[0], text_fields[1],
                    *coordinate_fields, text_fields[2], residue_id,
                )
            )
            if len(line) != 120:
                raise PreparationError("CHARMM EXT coordinate record is not 120 columns")
            stream.write(line + "\n")


def read_native_charmm_crd(path: Path) -> list[tuple[float, float, float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    header_index = next((index for index, line in enumerate(lines)
                         if line and not line.startswith("*")), None)
    if header_index is None:
        raise PreparationError("CHARMM EXT coordinate file has no header: %s" % path)
    header = lines[header_index]
    try:
        atom_count = int(header[:10])
    except ValueError as error:
        raise PreparationError("invalid CHARMM EXT atom count in %s: %s" % (path, error))
    if header[10:15] != "  EXT" or atom_count <= 0:
        raise PreparationError("invalid CHARMM EXT header in %s" % path)
    records = lines[header_index + 1:]
    if len(records) != atom_count:
        raise PreparationError("CHARMM EXT coordinate count differs from its header: %s" % path)
    coordinates = []
    for expected, line in enumerate(records, 1):
        if len(line) != 120:
            raise PreparationError("CHARMM EXT record is not 120 columns in %s" % path)
        try:
            atom_index = int(line[0:10])
            int(line[10:20])
            xyz = tuple(float(line[start:start + 20]) for start in (40, 60, 80))
            int(line[112:120])
        except ValueError as error:
            raise PreparationError("invalid CHARMM EXT coordinate in %s: %s" % (path, error))
        if atom_index != expected or any(not math.isfinite(value) for value in xyz):
            raise PreparationError("invalid CHARMM EXT atom order or coordinate in %s" % path)
        coordinates.append(xyz)
    return coordinates


def native_charmm_none_required(summary: dict) -> dict:
    charge = float(summary["raw_charge"])
    if integer_charge(charge, "native CHARMM source") != 0 or abs(charge) > POST_REPLACEMENT_TOLERANCE_E:
        raise PreparationError("native CHARMM archive is not formally neutral")
    return {
        "method": "none_required",
        "seed": None,
        "tolerance_e": FORMAL_CHARGE_TOLERANCE_E,
        "min_ion_separation_nm": None,
        "min_ion_solute_distance_nm": None,
        "neutralization_ion": None,
        "neutralization_ion_charge_e": None,
        "neutralization_ion_count": 0,
        "replaced_water_count": 0,
        "existing_ions": summary["existing_ions"],
        "water_count_before": summary["water_count"],
        "water_count_after": summary["water_count"],
        "selected_water_residue_indices": [],
        "selected_water_oxygen_positions_nm": [],
        "source_charge_e": charge,
        "source_formal_charge_e": 0,
        "pre_neutralization_charge_e": charge,
        "pre_neutralization_formal_charge_e": 0,
        "post_neutralization_charge_e": charge,
        "post_neutralization_formal_charge_e": 0,
    }


def export_native_charmm_archive(
    system: str, source_dir: Path, directory: Path,
) -> tuple[dict, dict, dict, dict]:
    specification = NATIVE_CHARMM_ARCHIVES[system]
    source_psf = source_dir / specification["psf"]
    source_pdb = source_dir / specification["pdb"]
    source_restart = source_dir / specification["restart"]
    topology_sources = tuple(source_dir / name for name in specification["topologies"])
    parameter_sources = tuple(source_dir / name for name in specification["parameters"])
    source_stream = source_dir / specification["stream"]
    for source in (
            source_psf, source_pdb, source_restart, source_stream,
            *topology_sources, *parameter_sources):
        if not source.is_file():
            raise PreparationError("native %s archive is missing %s" % (system, source))

    parsed = parse_surgical_charmm_psf(source_psf)
    summary = native_charmm_psf_summary(parsed)
    coordinates, box = read_tagged_genesis_restart(source_restart)
    if len(coordinates) != summary["atom_count"]:
        raise PreparationError("native %s restart/PSF atom counts differ" % system)

    shutil.copy2(source_psf, directory / "system.psf")
    shutil.copy2(source_pdb, directory / "system.pdb")
    write_native_charmm_crd(directory / "system.crd", parsed["atoms"], coordinates)
    write_namd_binary_coordinates(directory / "system.coor", coordinates)

    topology_names, flattened_parameter_names = flatten_native_charmm_stream(
        system, source_stream, topology_sources, directory,
    )
    parameter_names = []
    for source in parameter_sources:
        destination = directory / source.name
        if destination.exists():
            raise PreparationError("duplicate flattened CHARMM parameter name: %s" % source.name)
        shutil.copy2(source, destination)
        parameter_names.append(source.name)
    parameter_names.extend(flattened_parameter_names)
    topology_cards = tuple(directory / name for name in topology_names)
    write_native_charmm_xplor_psf(
        directory / "system.psf", directory / "system_xplor.psf", topology_cards,
    )
    hmr_validation = write_native_charmm_hmr_psf(
        directory / "system_xplor.psf", directory / "system_hmr_xplor.psf",
    )

    named = parse_surgical_charmm_psf(directory / "system_xplor.psf")
    defined_types = charmm_named_types(topology_cards)
    missing_types = sorted({atom["atomtype"] for atom in named["atoms"]} - defined_types)
    if missing_types:
        raise PreparationError("native %s PSF uses undefined named types: %s" % (
            system, ", ".join(missing_types),
        ))
    if parsed["lines"][parsed["atom_end"]:] != named["lines"][named["atom_end"]:]:
        raise PreparationError("numeric/XPLOR conversion changed native %s topology" % system)
    if len(parsed["atoms"]) != len(named["atoms"]):
        raise PreparationError("numeric/XPLOR conversion changed native %s atoms" % system)
    for original, converted in zip(parsed["atoms"], named["atoms"]):
        for field in ("index", "segid", "resid", "resname", "atomname", "charge", "mass"):
            if original[field] != converted[field]:
                raise PreparationError("numeric/XPLOR conversion changed %s atom metadata" % system)

    assets = {
        "GENESIS": {
            "format": "CHARMM",
            "topology": "system.psf",
            "coordinates": "system.crd",
            "topology_definitions": list(topology_names),
            "parameters": list(parameter_names),
        },
        "NAMD": {
            "format": "CHARMM",
            "topology": "system_xplor.psf",
            "topology_hmr": "system_hmr_xplor.psf",
            "coordinates": "system.coor",
            "coordinate_reference": "system.pdb",
            "parameters": list(parameter_names),
        },
    }
    context = {
        "coordinates": coordinates,
        "box": tuple(float(value) for value in box),
        "summary": summary,
        "hmr_validation": hmr_validation,
        "source_psf": source_psf,
        "source_pdb": source_pdb,
        "source_restart": source_restart,
        "source_stream": source_stream,
        "topology_sources": topology_sources,
        "parameter_sources": parameter_sources,
    }
    return assets, native_charmm_none_required(summary), {
        key: value for key, value in hmr_validation.items()
        if key in (
            "hydrogen_mass_dalton", "method", "particles_with_changed_mass",
            "water_particles_with_changed_mass", "total_mass_delta_dalton",
        )
    }, context


def validate_native_charmm_archive(
    system: str, directory: Path, assets: dict, context: dict,
) -> dict:
    import numpy as np

    specification = NATIVE_CHARMM_ARCHIVES[system]
    summary = context["summary"]
    expected_summary = {
        "atom_count": specification["atom_count"],
        "water_count": specification["water_count"],
        "existing_ions": specification["existing_ions"],
        "bonded_term_counts": specification["bonded_term_counts"],
        "exception_count": specification["exception_count"],
    }
    observed_summary = {key: summary[key] for key in expected_summary}
    if observed_summary != expected_summary:
        raise PreparationError("pinned native %s topology summary changed" % system)
    if summary["raw_charge"] != Decimal(0):
        raise PreparationError("pinned native %s PSF is not exactly neutral" % system)
    if summary["normal_mass"] != Decimal(specification["normal_total_mass_dalton"]):
        raise PreparationError("pinned native %s total mass changed" % system)
    if tuple(context["box"]) != tuple(specification["box_angstrom"]):
        raise PreparationError("pinned native %s restart box changed" % system)

    hmr = context["hmr_validation"]
    if (hmr["solute_hydrogen_count"] != specification["hmr_solute_hydrogens"]
            or hmr["donor_count"] != specification["hmr_donors"]
            or hmr["particles_with_changed_mass"]
            != specification["hmr_solute_hydrogens"] + specification["hmr_donors"]
            or hmr["normal_total_mass"] != summary["normal_mass"]
            or hmr["hmr_total_mass"] != summary["normal_mass"]):
        raise PreparationError("pinned native %s HMR validation changed" % system)

    crd_coordinates = np.asarray(read_native_charmm_crd(directory / "system.crd"), dtype=float)
    binary_coordinates = read_namd_binary_coordinates(directory / "system.coor")
    reference_coordinates = np.asarray(context["coordinates"], dtype=float)
    crd_delta = coordinate_delta(reference_coordinates, crd_coordinates)
    binary_delta = coordinate_delta(reference_coordinates, binary_coordinates)
    if (crd_delta > CHARMM_EXT_COORDINATE_TOLERANCE_ANGSTROM
            or binary_delta != 0.0):
        raise PreparationError(
            "native %s restart coordinates exceeded CHARMM EXT precision" % system
        )
    pdb_identities, _ = read_pdb_scaffold(directory / "system.pdb")
    if pdb_identities != summary["pdb_identities"]:
        raise PreparationError("native %s PDB scaffold changed PSF atom order/identity" % system)

    expected_hashes = specification["asset_sha256"]
    files = asset_paths(assets)
    if files != set(expected_hashes):
        raise PreparationError("native %s direct asset set differs from validated prototype" % system)
    observed_hashes = {name: sha256(directory / name) for name in sorted(files)}
    if observed_hashes != expected_hashes:
        mismatches = sorted(
            name for name in files if observed_hashes[name] != expected_hashes[name]
        )
        raise PreparationError("native %s prototype hashes changed: %s" % (
            system, ", ".join(mismatches),
        ))

    atom_count = summary["atom_count"]
    charge = float(summary["raw_charge"])
    mass = float(summary["normal_mass"])
    box = list(context["box"])
    canonical = {
        "atom_count": atom_count,
        "raw_charge_e": charge,
        "formal_charge_e": 0,
        "box_angstrom": box,
        "bonded_term_counts": summary["bonded_term_counts"],
        "exception_count": summary["exception_count"],
        "constraint_count": 0,
    }
    numerical = {
        "status": "not_run",
        "reference_variant": None,
        "energy_delta_kj_mol": None,
        "force_relative_rms": None,
        "cutoff_nm": None,
        "ewald_error_tolerance": None,
        "platform": None,
    }
    return {
        "status": "pass",
        "canonical": canonical,
        "assets": {
            "GENESIS": asset_record(atom_count, charge, box, crd_delta, mass, None),
            "NAMD": asset_record(atom_count, charge, box, binary_delta, mass, mass),
        },
        "representative_numerical": numerical,
    }


def asset_paths(assets: dict) -> set[str]:
    result = set()
    for record in assets.values():
        for key, value in record.items():
            if key == "format":
                continue
            if isinstance(value, str):
                result.add(value)
            elif isinstance(value, list):
                result.update(value)
    return result


def commit_directory(staging: Path, destination: Path, force: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        raise PreparationError("%s already exists (use --force to replace it)" % destination)
    backup = destination.with_name(".%s.backup-%d" % (destination.name, os.getpid()))
    if backup.exists():
        shutil.rmtree(backup)
    if destination.exists():
        os.replace(destination, backup)
    try:
        os.replace(staging, destination)
    except BaseException:
        if backup.exists() and not destination.exists():
            os.replace(backup, destination)
        raise
    if backup.exists():
        shutil.rmtree(backup)


def prepare_dhfr(
    variant: NativeVariant,
    output_root: Path,
    ion_parameters: Path,
    force: bool,
) -> Path:
    import openmm
    import parmed

    ion_parameter_sha256 = sha256(ion_parameters)
    if ion_parameter_sha256 != AMBERTOOLS26_TIP3P_ION_SHA256:
        raise PreparationError(
            "AmberTools26 TIP3P ion source SHA-256 changed: expected %s, got %s"
            % (AMBERTOOLS26_TIP3P_ION_SHA256, ion_parameter_sha256)
        )

    archive = DATA / "dhfr.tgz"
    with tempfile.TemporaryDirectory(prefix="genesis-benchmark-dhfr-") as temporary:
        source_dir = unpack_archive("dhfr", Path(temporary))
        source_topology = source_dir / "prmtop"
        source_coordinates = source_dir / "inpcrd"
        neutral, _, neutralization, preservation = mutate_native_dhfr(
            source_topology, source_coordinates, ion_parameters
        )
        hmr, hmr_metadata = copy_with_hmr(neutral)

        destination = output_root / variant.system / variant.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".%s.staging-" % variant.name, dir=destination.parent))
        try:
            assets, normal_system, hmr_system = export_native_dhfr(staging, neutral, hmr)
            validation = validate_native_dhfr(
                staging, neutral, normal_system, hmr_system, hmr_metadata, preservation
            )
            files = asset_paths(assets)
            missing = sorted(name for name in files if not (staging / name).is_file())
            if missing:
                raise PreparationError("missing exported assets: %s" % ", ".join(missing))
            checksums = {name: sha256(staging / name) for name in sorted(files)}
            manifest = {
                "schema": 1,
                "system": variant.system,
                "variant": variant.name,
                "family": variant.family,
                "forcefield": variant.forcefield,
                "water_model": variant.water_model,
                "solvent": "explicit",
                "atom_count": len(neutral.atoms),
                "box_angstrom": [float(value) for value in neutral.box[:3]],
                "assets": assets,
                "checksums": checksums,
                "hmr": hmr_metadata,
                "neutralization": neutralization,
                "provenance": {
                    "preparation_script": "prepare_variants.py",
                    "preparation_script_version": SCRIPT_VERSION,
                    "preparation_script_sha256": PREPARATION_SCRIPT_SHA256,
                    "source_transformations": {},
                    "source_archive": {"path": "data/dhfr.tgz", "sha256": sha256(archive)},
                    "source_atom_count": preservation["original_atom_count"],
                    "coordinate_source": {
                        "path": "data/dhfr.tgz::dhfr/inpcrd",
                        "sha256": sha256(source_coordinates),
                        "legacy_restart_used": False,
                    },
                    "forcefield_files": [
                        {"path": "data/dhfr.tgz::dhfr/prmtop", "sha256": sha256(source_topology)},
                        {"path": AMBERTOOLS26_TIP3P_ION_SOURCE,
                         "sha256": ion_parameter_sha256},
                    ],
                    "charmm_toppar_files": [],
                    "openmm_version": openmm.__version__,
                    "openmm_git_revision": openmm.version.git_revision,
                    "openmm_release": openmm.version.release,
                    "parmed_version": parmed.__version__,
                },
                "validation": validation,
            }
            (staging / "variant.json").write_text(stable_json(manifest), encoding="utf-8")
            commit_directory(staging, destination, force)
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
    return destination


def prepare_apoa1(
    variant: NativeVariant,
    output_root: Path,
    force: bool,
) -> Path:
    import openmm
    import parmed

    if (variant.family, variant.forcefield, variant.water_model) != (
            "CHARMM", "CHARMM27", "mTIP3P"):
        raise PreparationError("ApoA1 native identity must remain CHARMM27/mTIP3P")
    archive = DATA / "apoa1.tgz"
    with tempfile.TemporaryDirectory(prefix="genesis-benchmark-apoa1-") as temporary:
        source_dir = unpack_archive("apoa1", Path(temporary))
        source_topology = source_dir / "apoa1.psf"
        source_restart = source_dir / "apoa1.rst"
        source_rtf = source_dir / "top_all27_prot_lipid.rtf"
        source_prm = source_dir / "par_all27_prot_lipid.prm"
        for source_file in (source_topology, source_restart, source_rtf, source_prm):
            if not source_file.is_file():
                raise PreparationError("native ApoA1 archive is missing %s" % source_file.name)

        source, neutral, hmr, neutralization, hmr_metadata, preservation = (
            mutate_native_apoa1(source_topology, source_restart)
        )
        destination = output_root / variant.system / variant.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=".%s.staging-" % variant.name, dir=destination.parent,
        ))
        try:
            assets, parameters = export_native_apoa1(
                staging, neutral, hmr, source_rtf, source_prm,
            )
            validation, preservation_terms = validate_native_apoa1(
                staging, source, neutral, hmr_metadata, preservation, parameters,
            )
            expected_terms = {
                "bonds": {"source": 92118, "retained": 92076},
                "angles": {"source": 74136, "retained": 74122},
                "dihedrals": {"source": 74130, "retained": 74130},
                "impropers": {"source": 1402, "retained": 1402},
                "donors": {"source": 758, "retained": 758},
                "acceptors": {"source": 22098, "retained": 22084},
                "cmaps": {"source": 388, "retained": 388},
            }
            if preservation_terms != expected_terms:
                raise PreparationError("pinned ApoA1 native topology counts changed")
            files = asset_paths(assets)
            missing = sorted(name for name in files if not (staging / name).is_file())
            if missing:
                raise PreparationError("missing exported ApoA1 assets: %s" % ", ".join(missing))
            extra = sorted(
                path.name for path in staging.iterdir()
                if path.is_file() and path.name not in files
            )
            if extra:
                raise PreparationError("unadvertised ApoA1 assets were generated: %s" % (
                    ", ".join(extra),
                ))
            checksums = {name: sha256(staging / name) for name in sorted(files)}
            manifest = {
                "schema": 1,
                "system": variant.system,
                "variant": variant.name,
                "family": variant.family,
                "forcefield": variant.forcefield,
                "water_model": variant.water_model,
                "solvent": "explicit",
                "atom_count": len(neutral.atoms),
                "box_angstrom": [float(value) for value in neutral.box[:3]],
                "assets": assets,
                "checksums": checksums,
                "hmr": hmr_metadata,
                "neutralization": neutralization,
                "provenance": {
                    "preparation_script": "prepare_variants.py",
                    "preparation_script_version": SCRIPT_VERSION,
                    "preparation_script_sha256": PREPARATION_SCRIPT_SHA256,
                    "source_transformations": {},
                    "source_archive": {"path": "data/apoa1.tgz", "sha256": sha256(archive)},
                    "source_atom_count": preservation["source_atom_count"],
                    "coordinate_source": {
                        "path": "data/apoa1.tgz::apoa1/apoa1.rst",
                        "sha256": sha256(source_restart),
                        "legacy_restart_used": True,
                    },
                    "forcefield_files": [
                        {"path": "data/apoa1.tgz::apoa1/apoa1.psf",
                         "sha256": sha256(source_topology)},
                    ],
                    "charmm_toppar_files": [
                        {"path": "data/apoa1.tgz::apoa1/top_all27_prot_lipid.rtf",
                         "sha256": sha256(source_rtf)},
                        {"path": "data/apoa1.tgz::apoa1/par_all27_prot_lipid.prm",
                         "sha256": sha256(source_prm)},
                    ],
                    "openmm_version": openmm.__version__,
                    "openmm_git_revision": openmm.version.git_revision,
                    "openmm_release": openmm.version.release,
                    "parmed_version": parmed.__version__,
                },
                "validation": validation,
            }
            (staging / "variant.json").write_text(stable_json(manifest), encoding="utf-8")
            commit_directory(staging, destination, force)
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
    return destination


def prepare_native_amber_archive(
    variant: NativeVariant,
    output_root: Path,
    force: bool,
) -> Path:
    """Prepare a neutral source-preserving AMBER model and static HMR."""
    import openmm
    import parmed

    try:
        specification = NATIVE_AMBER_ARCHIVES[variant.system]
    except KeyError as error:
        raise PreparationError("no native Amber archive specification for %s" % variant.system) from error
    if variant.family != "AMBER" or variant.water_model != "TIP3P":
        raise PreparationError("native Amber archive identity changed for %s" % variant.system)

    archive = DATA / (variant.system + ".tgz")
    with tempfile.TemporaryDirectory(
            prefix="genesis-benchmark-%s-" % variant.system) as temporary:
        source_dir = unpack_archive(variant.system, Path(temporary))
        source_topology = source_dir / specification["topology"]
        source_coordinates = source_dir / specification["coordinates"]
        for source_file in (source_topology, source_coordinates):
            if not source_file.is_file():
                raise PreparationError(
                    "native %s archive is missing %s" % (variant.system, source_file.name)
                )

        destination = output_root / variant.system / variant.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=".%s.staging-" % variant.name, dir=destination.parent,
        ))
        try:
            shutil.copyfile(source_topology, staging / "system.prmtop")
            shutil.copyfile(source_coordinates, staging / "system.inpcrd")
            dihedral_normalization = canonicalize_native_amber_multiterm_dihedrals(
                staging / "system.prmtop",
            )
            box_reconciliation = reconcile_native_amber_box(
                staging / "system.prmtop", staging / "system.inpcrd",
            )
            normal = parmed.load_file(
                str(staging / "system.prmtop"), xyz=str(staging / "system.inpcrd")
            )
            if len(normal.atoms) != specification["atom_count"]:
                raise PreparationError(
                    "%s native atom count changed: expected %d, got %d" % (
                        variant.system, specification["atom_count"], len(normal.atoms),
                    )
                )
            if normal.box is None or any(
                    abs(float(observed) - expected) > 1.0e-7
                    for observed, expected in zip(
                        normal.box[:3], specification["box_angstrom"],
                    )):
                raise PreparationError("%s authoritative restart box changed" % variant.system)
            neutralization = native_amber_none_required(normal)
            hmr, hmr_metadata = copy_with_hmr(normal)

            assets, normal_system, runtime_hmr_system = export_native_amber_archive(
                staging, normal, hmr, specification["gromacs_equivalence"],
            )
            validation = validate_native_amber_archive(
                staging, source_topology, source_coordinates, normal,
                normal_system, runtime_hmr_system, hmr_metadata, assets,
                box_reconciliation, dihedral_normalization, variant.name,
            )
            files = asset_paths(assets)
            missing = sorted(name for name in files if not (staging / name).is_file())
            if missing:
                raise PreparationError(
                    "missing exported %s assets: %s" % (
                        variant.system, ", ".join(missing),
                    )
                )
            extra = sorted(
                path.name for path in staging.iterdir()
                if path.is_file() and path.name not in files
            )
            if extra:
                raise PreparationError(
                    "unadvertised %s assets were generated: %s" % (
                        variant.system, ", ".join(extra),
                    )
                )
            checksums = {name: sha256(staging / name) for name in sorted(files)}
            member_prefix = "data/%s.tgz::%s/" % (variant.system, variant.system)
            manifest = {
                "schema": 1,
                "system": variant.system,
                "variant": variant.name,
                "family": variant.family,
                "forcefield": variant.forcefield,
                "water_model": variant.water_model,
                "solvent": "explicit",
                "atom_count": len(normal.atoms),
                "box_angstrom": [float(value) for value in normal.box[:3]],
                "assets": assets,
                "checksums": checksums,
                "hmr": hmr_metadata,
                "neutralization": neutralization,
                "provenance": {
                    "preparation_script": "prepare_variants.py",
                    "preparation_script_version": SCRIPT_VERSION,
                    "preparation_script_sha256": PREPARATION_SCRIPT_SHA256,
                    "source_transformations": native_amber_source_transformations(
                        box_reconciliation, dihedral_normalization,
                    ),
                    "source_archive": {
                        "path": "data/%s.tgz" % variant.system,
                        "sha256": sha256(archive),
                    },
                    "source_atom_count": len(normal.atoms),
                    "coordinate_source": {
                        "path": member_prefix + specification["coordinates"],
                        "sha256": sha256(source_coordinates),
                        "legacy_restart_used": False,
                    },
                    "forcefield_files": [{
                        "path": member_prefix + specification["topology"],
                        "sha256": sha256(source_topology),
                    }],
                    "charmm_toppar_files": [],
                    "openmm_version": openmm.__version__,
                    "openmm_git_revision": openmm.version.git_revision,
                    "openmm_release": openmm.version.release,
                    "parmed_version": parmed.__version__,
                },
                "validation": validation,
            }
            (staging / "variant.json").write_text(stable_json(manifest), encoding="utf-8")
            commit_directory(staging, destination, force)
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
    return destination


def prepare_native_charmm_archive(
    variant: NativeVariant,
    output_root: Path,
    force: bool,
) -> Path:
    import openmm
    import parmed

    if variant.system not in NATIVE_CHARMM_ARCHIVES:
        raise PreparationError("no source-preserving CHARMM archive builder for %s" % variant.system)
    specification = NATIVE_CHARMM_ARCHIVES[variant.system]
    if (variant.family, variant.forcefield, variant.water_model) != (
            "CHARMM", specification["forcefield"], "mTIP3P"):
        raise PreparationError("%s native CHARMM identity changed" % variant.system)
    archive = DATA / (variant.system + ".tgz")
    with tempfile.TemporaryDirectory(
            prefix="genesis-benchmark-%s-" % variant.system) as temporary:
        source_dir = unpack_archive(variant.system, Path(temporary))
        destination = output_root / variant.system / variant.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=".%s.staging-" % variant.name, dir=destination.parent,
        ))
        try:
            assets, neutralization, hmr, context = export_native_charmm_archive(
                variant.system, source_dir, staging,
            )
            validation = validate_native_charmm_archive(
                variant.system, staging, assets, context,
            )
            files = asset_paths(assets)
            missing = sorted(name for name in files if not (staging / name).is_file())
            if missing:
                raise PreparationError("missing exported %s assets: %s" % (
                    variant.system, ", ".join(missing),
                ))
            extra = sorted(
                path.name for path in staging.iterdir()
                if path.is_file() and path.name not in files
            )
            if extra:
                raise PreparationError("unadvertised %s assets were generated: %s" % (
                    variant.system, ", ".join(extra),
                ))
            checksums = {name: sha256(staging / name) for name in sorted(files)}
            member = lambda name: "data/%s.tgz::%s/%s" % (
                variant.system, variant.system, name,
            )
            source_restart = context["source_restart"]
            charmm_sources = (
                tuple(context["topology_sources"])
                + tuple(context["parameter_sources"])
                + (context["source_stream"],)
            )
            manifest = {
                "schema": 1,
                "system": variant.system,
                "variant": variant.name,
                "family": variant.family,
                "forcefield": variant.forcefield,
                "water_model": variant.water_model,
                "solvent": "explicit",
                "atom_count": context["summary"]["atom_count"],
                "box_angstrom": list(context["box"]),
                "assets": assets,
                "checksums": checksums,
                "hmr": hmr,
                "neutralization": neutralization,
                "provenance": {
                    "preparation_script": "prepare_variants.py",
                    "preparation_script_version": SCRIPT_VERSION,
                    "preparation_script_sha256": PREPARATION_SCRIPT_SHA256,
                    "source_transformations": {},
                    "source_archive": {
                        "path": "data/%s.tgz" % variant.system,
                        "sha256": sha256(archive),
                    },
                    "source_atom_count": context["summary"]["atom_count"],
                    "coordinate_source": {
                        "path": member(specification["restart"]),
                        "sha256": sha256(source_restart),
                        "legacy_restart_used": True,
                    },
                    "forcefield_files": [{
                        "path": member(specification["psf"]),
                        "sha256": sha256(context["source_psf"]),
                    }],
                    "charmm_toppar_files": [
                        {
                            "path": member(str(source.relative_to(source_dir))),
                            "sha256": sha256(source),
                        }
                        for source in charmm_sources
                    ],
                    "openmm_version": openmm.__version__,
                    "openmm_git_revision": openmm.version.git_revision,
                    "openmm_release": openmm.version.release,
                    "parmed_version": parmed.__version__,
                },
                "validation": validation,
            }
            (staging / "variant.json").write_text(stable_json(manifest), encoding="utf-8")
            commit_directory(staging, destination, force)
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
    return destination


def select_variants(args: argparse.Namespace) -> list[NativeVariant]:
    systems = parse_csv(args.systems)
    selected = [
        variant for variant in VARIANTS
        if (systems is None or variant.system in systems)
    ]
    unknown = sorted((systems or set()) - {variant.system for variant in VARIANTS})
    if unknown:
        raise PreparationError("unknown system: %s" % ", ".join(unknown))
    if not selected:
        raise PreparationError("no native variants match the requested filters")
    return selected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--systems", action="append",
        help="native system name(s), repeatable or comma-separated",
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--amber-tip3p-ions", type=Path, default=DEFAULT_AMBER_TIP3P_IONS)
    parser.add_argument("--force", action="store_true", help="atomically replace an existing prepared variant")
    parser.add_argument("--dry-run", action="store_true", help="print the selected action without loading data")
    parser.add_argument(
        "--list", action="store_true",
        help="list the one canonical native model per system",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        selected = select_variants(args)
        if args.list:
            for variant in selected:
                suffix = "" if not variant.reason else " - " + variant.reason
                print("%-10s %-18s %-8s %-8s %s%s" % (
                    variant.system, variant.forcefield, variant.water_model,
                    variant.family, variant.status, suffix,
                ))
            return 0

        actionable = [variant for variant in selected if variant.status == "ready"]
        deferred = [variant for variant in selected if variant.status == "deferred"]
        archive = [variant for variant in selected if variant.status == "archive"]
        if deferred:
            details = "; ".join("%s: %s" % (v.system, v.reason) for v in deferred)
            raise PreparationError("requested native transformation is unfinished: %s" % details)
        for variant in archive:
            print("%s: no topology edit required; canonical archive fallback remains active"
                  % variant.system)
        for variant in actionable:
            if variant.system in NATIVE_CHARMM_ARCHIVES:
                action = "export restart-derived native CHARMM direct assets"
            elif variant.system in NATIVE_AMBER_ARCHIVES:
                action = "export source-preserving native AMBER direct/HMR assets"
            else:
                action = "replace native waters with native-family counterions"
            print("%s: %s" % (variant.system, action))
        if args.dry_run or not actionable:
            return 0

        try:
            import openmm  # noqa: F401
            import parmed  # noqa: F401
        except ImportError as error:
            raise PreparationError(
                "OpenMM and ParmEd are required; run with the OpenMM conda environment: %s" % error
            )
        toolchain = (
            openmm.__version__, openmm.version.git_revision,
            openmm.version.release, parmed.__version__,
        )
        expected_toolchain = (
            EXPECTED_OPENMM_VERSION, EXPECTED_OPENMM_GIT_REVISION,
            EXPECTED_OPENMM_RELEASE, EXPECTED_PARMED_VERSION,
        )
        if toolchain != expected_toolchain:
            raise PreparationError(
                "prepared assets require the pinned OpenMM/ParmEd toolchain %r; got %r"
                % (expected_toolchain, toolchain)
            )
        for variant in actionable:
            if variant.system == "dhfr":
                destination = prepare_dhfr(
                    variant, args.output_root.resolve(),
                    args.amber_tip3p_ions.resolve(), args.force,
                )
            elif variant.system == "apoa1":
                destination = prepare_apoa1(
                    variant, args.output_root.resolve(), args.force,
                )
            elif variant.system in NATIVE_AMBER_ARCHIVES:
                destination = prepare_native_amber_archive(
                    variant, args.output_root.resolve(), args.force,
                )
            elif variant.system in NATIVE_CHARMM_ARCHIVES:
                destination = prepare_native_charmm_archive(
                    variant, args.output_root.resolve(), args.force,
                )
            else:
                raise PreparationError("no enabled builder for %s" % variant.system)
            print("wrote %s" % destination)
        return 0
    except PreparationError as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    sys.exit(main())
