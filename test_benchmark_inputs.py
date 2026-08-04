import os
import re
import tarfile
import unittest

import generate_inputs
import run_benchmark


RETIRED_KEYS = {
    "cell_size_autotune",
    "group_tp",
    "nbupdate_autotune",
    "nbupdate_period",
    "nonbond_kernel",
    "water_model",
}

AMBER_TOPOLOGIES = {
    "ake": ("ake.tgz", "ake/ake.top"),
    "cellulose": ("cellulose.tgz", "cellulose/prmtop"),
    "dhfr": ("dhfr.tgz", "dhfr/prmtop"),
    "factorix": ("factorix.tgz", "factorix/FactorIX.prmtop"),
    "stmv": ("stmv.tgz", "stmv/prmtop"),
}

PSF_TOPOLOGIES = {
    "apoa1": ("apoa1.tgz", "apoa1/apoa1.psf", 92224, 59232),
    "dppc": ("dppc.tgz", "dppc/dppc.psf", 36126, 23004),
    "uun": ("uun.tgz", "uun/uun.psf", 216726, 137196),
}


def amber_flag_values(topology, flag):
    match = re.search(
        r"^%FLAG "
        + re.escape(flag)
        + r"\s*\n%FORMAT\((\d+)([A-Za-z])(\d+)(?:\.\d+)?\)\s*\n"
        + r"(.*?)(?=^%FLAG |\Z)",
        topology,
        re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise AssertionError("missing AMBER topology flag %s" % flag)
    width = int(match.group(3))
    values = []
    for line in match.group(4).splitlines():
        values.extend(
            line[index : index + width].strip()
            for index in range(0, len(line), width)
            if line[index : index + width].strip()
        )
    return values


def psf_atoms(topology):
    lines = topology.splitlines()
    header_index = next(index for index, line in enumerate(lines) if "!NATOM" in line)
    atom_count = int(lines[header_index].split()[0])
    fields = [line.split() for line in lines[header_index + 1 : header_index + 1 + atom_count]]
    return [(atom[4], float(atom[7])) for atom in fields]


def gromacs_atoms(topology):
    atoms = []
    in_atoms = False
    for raw_line in topology.splitlines():
        line = raw_line.split(";", 1)[0].strip()
        section = re.match(r"^\[\s*([^]]+)\s*\]$", line)
        if section:
            in_atoms = section.group(1).strip().lower() == "atoms"
            continue
        if line.startswith("#"):
            continue
        if in_atoms and line:
            fields = line.split()
            atoms.append((fields[4], float(fields[7])))
    return atoms


def assignments(text):
    result = {}
    section = None
    for line in text.splitlines():
        section_match = run_benchmark.SECTION_RE.match(line)
        if section_match:
            section = section_match.group(1).upper()
            continue
        key_match = run_benchmark.KEYLINE_RE.match(line)
        if key_match:
            result[(section, key_match.group(1).lower())] = line.split("=", 1)[1].strip()
    return result


class BenchmarkInputTests(unittest.TestCase):
    def test_generated_matrix_is_complete_and_current(self):
        self.assertEqual(list(generate_inputs.SYSTEMS), run_benchmark.ALL_SYSTEMS)
        expected_names = {
            "%s_%s_%s.inp" % (system, ensemble, time_step)
            for system in run_benchmark.ALL_SYSTEMS
            for ensemble in run_benchmark.ALL_ENSEMBLES
            for time_step in run_benchmark.ALL_TIME_STEPS
        }
        actual_names = {
            name for name in os.listdir(generate_inputs.OUTPUT_DIR) if name.endswith(".inp")
        }
        self.assertEqual(len(expected_names), 54)
        self.assertEqual(actual_names, expected_names)

        for name in expected_names:
            system, ensemble, time_step = name[:-4].rsplit("_", 2)
            path = os.path.join(generate_inputs.OUTPUT_DIR, name)
            with open(path) as input_file:
                text = input_file.read()
            self.assertEqual(text, generate_inputs.make_input(system, ensemble, time_step))

    def test_assignments_are_aligned_and_retired_keys_are_absent(self):
        for name in os.listdir(generate_inputs.OUTPUT_DIR):
            if not name.endswith(".inp"):
                continue
            with self.subTest(input=name):
                with open(os.path.join(generate_inputs.OUTPUT_DIR, name)) as input_file:
                    text = input_file.read()
                assignment_lines = [
                    line for line in text.splitlines() if run_benchmark.KEYLINE_RE.match(line)
                ]
                equals_columns = {line.index("=") for line in assignment_lines}
                keys = {
                    run_benchmark.KEYLINE_RE.match(line).group(1).lower()
                    for line in assignment_lines
                }
                self.assertEqual(len(equals_columns), 1)
                self.assertTrue(RETIRED_KEYS.isdisjoint(keys))

    def test_every_input_uses_requested_constraint_policy(self):
        for name in os.listdir(generate_inputs.OUTPUT_DIR):
            if not name.endswith(".inp"):
                continue
            with self.subTest(input=name):
                with open(os.path.join(generate_inputs.OUTPUT_DIR, name)) as input_file:
                    values = assignments(input_file.read())
                self.assertEqual(values[("CONSTRAINTS", "rigid_bond")], "YES")
                self.assertEqual(values[("CONSTRAINTS", "cons_scheme")], "MSHAKE")
                self.assertEqual(values[("CONSTRAINTS", "iter_solute")], "3")
                self.assertEqual(values[("CONSTRAINTS", "iter_water")], "3")

    def test_hmr_policy_uses_runtime_scaling_only_at_4fs(self):
        for system in run_benchmark.ALL_SYSTEMS:
            self.assertNotIn("hmr_topology", generate_inputs.SYSTEMS[system])
            for time_step in run_benchmark.ALL_TIME_STEPS:
                values = assignments(generate_inputs.make_input(system, "nve", time_step))
                has_runtime_hmr = ("DYNAMICS", "hydrogen_mr") in values
                has_mass_bound = ("CONSTRAINTS", "hydrogen_mass_upper_bound") in values
                if time_step == "4fs":
                    self.assertTrue(has_runtime_hmr)
                    self.assertTrue(has_mass_bound)
                    self.assertEqual(values[("DYNAMICS", "hmr_target")], "all")
                    self.assertEqual(values[("DYNAMICS", "hmr_ratio")], "3.0")
                else:
                    self.assertFalse(has_runtime_hmr)
                    self.assertFalse(has_mass_bound)

        factorix_inputs = generate_inputs.SYSTEMS["factorix"]["input"]
        self.assertFalse(
            any(line.split("=", 1)[0].strip() == "rstfile" for line in factorix_inputs)
        )
        factorix_npt = assignments(generate_inputs.make_input("factorix", "npt", "2fs"))
        self.assertEqual(factorix_npt[("BOUNDARY", "box_size_x")], "142.0855468")
        self.assertEqual(factorix_npt[("BOUNDARY", "box_size_y")], "83.3368905")
        self.assertEqual(factorix_npt[("BOUNDARY", "box_size_z")], "78.6783548")

    def test_all_amber_topologies_use_normal_hydrogen_masses(self):
        expected = {
            "factorix": (90906, 59533, 554212.334049),
            "stmv": (1067095, 677955, 6695311.430000),
        }
        for system, (archive_name, member_name) in AMBER_TOPOLOGIES.items():
            with self.subTest(system=system):
                archive_path = os.path.join(run_benchmark.DATA, archive_name)
                with tarfile.open(archive_path, "r:gz") as archive:
                    topology = archive.extractfile(member_name).read().decode("ascii")
                names = amber_flag_values(topology, "ATOM_NAME")
                masses = [float(value) for value in amber_flag_values(topology, "MASS")]
                self.assertEqual(len(names), len(masses))
                hydrogen_masses = [
                    mass for name, mass in zip(names, masses) if name.upper().startswith("H")
                ]
                self.assertTrue(hydrogen_masses)
                self.assertTrue(all(abs(mass - 1.008) < 1.0e-8 for mass in hydrogen_masses))
                if system in expected:
                    atom_count, hydrogen_count, total_mass = expected[system]
                    self.assertEqual(len(masses), atom_count)
                    self.assertEqual(len(hydrogen_masses), hydrogen_count)
                    self.assertAlmostEqual(sum(masses), total_mass, places=5)

    def test_all_psf_and_gromacs_topologies_use_normal_hydrogen_masses(self):
        for system, (archive_name, member_name, atom_count, hydrogen_count) in PSF_TOPOLOGIES.items():
            with self.subTest(system=system):
                archive_path = os.path.join(run_benchmark.DATA, archive_name)
                with tarfile.open(archive_path, "r:gz") as archive:
                    topology = archive.extractfile(member_name).read().decode("ascii")
                atoms = psf_atoms(topology)
                hydrogen_masses = [mass for name, mass in atoms if name.upper().startswith("H")]
                self.assertEqual(len(atoms), atom_count)
                self.assertEqual(len(hydrogen_masses), hydrogen_count)
                self.assertTrue(all(abs(mass - 1.008) < 1.0e-8 for mass in hydrogen_masses))

        with tarfile.open(os.path.join(run_benchmark.DATA, "bpti.tgz"), "r:gz") as archive:
            topology = archive.extractfile("bpti/bpti.top").read().decode("ascii")
            water_topology = archive.extractfile("bpti/amber03.ff/tip3p.itp").read().decode("ascii")
        protein_atoms = gromacs_atoms(topology)
        water_atoms = gromacs_atoms(water_topology)
        hydrogen_masses = [
            mass
            for name, mass in protein_atoms + water_atoms
            if name.upper().startswith("H")
        ]
        self.assertEqual(len(protein_atoms) + 8938 * len(water_atoms) + 6, 27712)
        self.assertEqual(
            sum(name.upper().startswith("H") for name, _ in protein_atoms)
            + 8938 * sum(name.upper().startswith("H") for name, _ in water_atoms),
            18314,
        )
        self.assertTrue(all(abs(mass - 1.008) < 1.0e-8 for mass in hydrogen_masses))

    def test_dhfr_archive_is_the_23558_atom_jac_system(self):
        archive_path = os.path.join(run_benchmark.DATA, "dhfr.tgz")
        with tarfile.open(archive_path, "r:gz") as archive:
            self.assertEqual(
                {member.name for member in archive.getmembers()},
                {"dhfr", "dhfr/inpcrd", "dhfr/prmtop"},
            )
            topology = archive.extractfile("dhfr/prmtop").read().decode("ascii")

        pointer_match = re.search(
            r"%FLAG POINTERS\s+%FORMAT\([^\n]+\)\s+([0-9]+)", topology
        )
        self.assertIsNotNone(pointer_match)
        self.assertEqual(int(pointer_match.group(1)), 23558)

    def test_driver_generated_inputs_remain_aligned_and_current(self):
        base = generate_inputs.make_input("dhfr", "npt", "4fs")
        tuned = {"kernel": dict.fromkeys(run_benchmark.KERNEL_KEYS, 128)}
        tuned["pme_schedule"] = "OVERLAP_BONDED"
        tune_text = run_benchmark.build_tune_input(base, {"kernel"}, 1000, 100, 5)
        pinned_text = run_benchmark.build_pinned_input(
            base, tuned, {"kernel"}, 1000, 100, 5
        )
        for text in (tune_text, pinned_text):
            lines = [line for line in text.splitlines() if run_benchmark.KEYLINE_RE.match(line)]
            self.assertEqual(len({line.index("=") for line in lines}), 1)
            keys = {run_benchmark.KEYLINE_RE.match(line).group(1).lower() for line in lines}
            self.assertTrue(RETIRED_KEYS.isdisjoint(keys))
        self.assertEqual(
            assignments(pinned_text)[("GPU", "pme_schedule")], "OVERLAP_BONDED"
        )

    def test_current_kernel_autotune_report_is_fully_parsed(self):
        report = "\n".join(
            "%s = %d" % (key, 64 + index * 32)
            for index, key in enumerate(run_benchmark.KERNEL_KEYS)
        ) + "\npme_schedule = OVERLAP_NONBONDED"
        tuned = run_benchmark.parse_autotune_report(report)
        self.assertEqual(set(tuned["kernel"]), set(run_benchmark.KERNEL_KEYS))
        self.assertEqual(tuned["pme_schedule"], "OVERLAP_NONBONDED")
        self.assertEqual(run_benchmark.validate_tuned_values(tuned, {"kernel"}), [])


if __name__ == "__main__":
    unittest.main()
