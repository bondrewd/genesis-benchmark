#!/usr/bin/env python3
"""Run GENESIS spdyn GPU benchmarks.

For each selected system, ensemble, and time step, the driver writes a run
input, warms up, measures production runs, and records raw logs plus a
row-per-run CSV.
"""

import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import shlex
import shutil
import signal
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time

import generate_inputs as input_generator


HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
INPUTS = os.path.join(HERE, "inputs", "GENESIS")
RESULTS = os.path.join(HERE, "results")
SPDYN = os.path.normpath(os.path.join(
    HERE, "..", "genesis-mkl-private-gpu", "src", "spdyn_singlempi", "spdyn",
))

ALL_SYSTEMS = list(input_generator.SYSTEM_ORDER)
ALL_ENSEMBLES = list(input_generator.ENSEMBLE_ORDER)
ARCHIVE_SENTINEL = ".archive_sha256"
PERFORMANCE_RE = re.compile(r"\[PERFORMANCE\]\s*(?:performance:\s*)?([0-9.]+)\s*ns/day")
GPU_ROUTE_PATTERNS = (
    ("real-space backend", re.compile(r"real_space_backend\s*=\s*GPU", re.IGNORECASE)),
    ("reciprocal-space backend", re.compile(r"reciprocal_space_backend\s*=\s*GPU", re.IGNORECASE)),
    ("pair-list route", re.compile(r"Pairlist\s*=\s*GPU", re.IGNORECASE)),
    ("nonbonded route", re.compile(r"Nonbond\s*=\s*GPU", re.IGNORECASE)),
)
NUM_ATOMS_PATTERNS = [
    re.compile(r"^\s*num_atoms\s*=\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
    re.compile(r"^\s*natom\s*=\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
    re.compile(r"^\s*number\s+of\s+atoms\s*[:=]\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
]
SECTION_RE = re.compile(r"^\s*\[([A-Za-z0-9_]+)\]\s*$")
KEYLINE_RE = re.compile(r"^\s*([A-Za-z_]\w*)\s*=")

BENCHMARK_LOG_FILE = None
BENCHMARK_LOG_BUFFER = []

# These are the only base-input controls a benchmark invocation may override.
# Every other option must match a fresh canonical render byte-for-value after
# parsing, so topology paths, algorithms, HMR, PME, constraints, and switching
# semantics cannot drift unnoticed.
MUTABLE_BASE_OPTIONS = {
    "input_energy_switchdist",
    "input_energy_cutoffdist",
    "input_energy_pairlistdist",
    "input_dynamics_iseed",
    "input_dynamics_nsteps",
    "input_dynamics_eneout_period",
    "input_dynamics_thermostat_period",
    "input_dynamics_barostat_period",
    "input_dynamics_baroscale_period",
    "input_ensemble_temperature",
    "input_ensemble_tau_t",
    "input_ensemble_pressure",
    "input_ensemble_tau_p",
}


def split_sections(text):
    """Split a GENESIS control file into section blocks.

    The return value is a list of two-item lists: section name, then the raw
    lines belonging to that section. Text before the first section has a
    section name of None.
    """
    blocks = []
    section_name = None
    section_lines = []
    for raw_line in text.splitlines():
        section_match = SECTION_RE.match(raw_line)
        if section_match:
            blocks.append([section_name, section_lines])
            section_name = section_match.group(1)
            section_lines = [raw_line]
        else:
            section_lines.append(raw_line)
    blocks.append([section_name, section_lines])
    return blocks


def join_sections(blocks):
    """Join section blocks into a GENESIS control-file string."""
    output_lines = []
    for section_name, section_lines in blocks:
        del section_name
        output_lines.extend(section_lines)
    return align_assignments("\n".join(output_lines))


def align_assignments(text):
    """Return control-file text with every parameter assignment aligned."""
    lines = text.splitlines()
    keys = []
    for raw_line in lines:
        key_match = KEYLINE_RE.match(raw_line)
        if key_match:
            keys.append(key_match.group(1))
    if not keys:
        return text.rstrip("\n") + "\n"

    width = max(len(key) for key in keys)
    aligned_lines = []
    for raw_line in lines:
        key_match = KEYLINE_RE.match(raw_line)
        if not key_match:
            aligned_lines.append(raw_line)
            continue
        key = key_match.group(1)
        value = raw_line.split("=", 1)[1].strip()
        aligned_lines.append("%s = %s" % (key.ljust(width), value))
    return "\n".join(aligned_lines).rstrip("\n") + "\n"


def find_section(blocks, section):
    """Return the section block whose name matches section, ignoring case."""
    for block in blocks:
        block_name = block[0]
        if block_name and block_name.upper() == section.upper():
            return block
    return None


def set_parameter(blocks, section, key, value):
    """Set a GENESIS key in a section, creating the section or key if needed."""
    block = find_section(blocks, section)
    formatted_value = str(value)
    if block is None:
        blocks.append([section, ["[%s]" % section, "%-16s = %s" % (key, formatted_value)]])
        return

    section_lines = block[1]
    previous_line_continues = False
    for line_index in range(1, len(section_lines)):
        raw_line = section_lines[line_index]
        stripped_line = raw_line.strip()
        continued_line = previous_line_continues
        previous_line_continues = stripped_line.endswith("\\")
        if continued_line or stripped_line.startswith("#"):
            continue
        key_match = KEYLINE_RE.match(raw_line)
        if key_match and key_match.group(1).lower() == key.lower():
            section_lines[line_index] = "%-16s = %s" % (key, formatted_value)
            return
    section_lines.insert(1, "%-16s = %s" % (key, formatted_value))


def has_parameter(blocks, section, key):
    """Return True when a section contains a key assignment."""
    block = find_section(blocks, section)
    if block is None:
        return False

    previous_line_continues = False
    for raw_line in block[1][1:]:
        stripped_line = raw_line.strip()
        continued_line = previous_line_continues
        previous_line_continues = stripped_line.endswith("\\")
        if continued_line or stripped_line.startswith("#"):
            continue
        key_match = KEYLINE_RE.match(raw_line)
        if key_match and key_match.group(1).lower() == key.lower():
            return True
    return False


def get_parameter_value(text, section, key):
    """Return an uncommented GENESIS key value from a section."""
    block = find_section(split_sections(text), section)
    if block is None:
        return None

    previous_line_continues = False
    for raw_line in block[1][1:]:
        stripped_line = raw_line.strip()
        continued_line = previous_line_continues
        previous_line_continues = stripped_line.endswith("\\")
        if continued_line or stripped_line.startswith("#"):
            continue
        key_match = KEYLINE_RE.match(raw_line)
        if key_match and key_match.group(1).lower() == key.lower():
            return raw_line.split("=", 1)[1].split("#", 1)[0].strip()
    return None


def set_run_window(blocks, num_steps, eneout_period):
    """Set the number of MD steps and the energy-output period."""
    set_parameter(blocks, "DYNAMICS", "nsteps", num_steps)
    set_parameter(blocks, "DYNAMICS", "eneout_period", eneout_period)


def float_parameter(blocks, section, key):
    """Return a finite numeric parameter from section blocks."""
    value = get_parameter_value(join_sections(blocks), section, key)
    if value is None:
        raise RuntimeError("input is missing [%s] %s" % (section, key))
    try:
        result = float(value)
    except ValueError:
        raise RuntimeError("[%s] %s is not numeric: %r" % (section, key, value))
    if not math.isfinite(result):
        raise RuntimeError("[%s] %s is not finite" % (section, key))
    return result


def integer_parameter(blocks, section, key):
    """Return an integer-valued control parameter from section blocks."""
    value = float_parameter(blocks, section, key)
    integer = int(round(value))
    if abs(value - integer) > 1.0e-12:
        raise RuntimeError("[%s] %s must be an integer" % (section, key))
    return integer


def apply_protocol_overrides(blocks, variant, args):
    """Apply requested GENESIS protocol overrides and validate the result."""
    set_run_window(blocks, args.num_steps, args.eneout_period)
    has_thermostat = has_parameter(blocks, "DYNAMICS", "thermostat_period")
    has_barostat = has_parameter(blocks, "DYNAMICS", "barostat_period")
    if args.thermostat_period is not None and has_thermostat:
        set_parameter(blocks, "DYNAMICS", "thermostat_period", args.thermostat_period)
    if args.barostat_period is not None and has_barostat:
        set_parameter(blocks, "DYNAMICS", "barostat_period", args.barostat_period)
    if args.baroscale_period is not None and has_barostat:
        set_parameter(blocks, "DYNAMICS", "baroscale_period", args.baroscale_period)
    if args.seed is not None:
        set_parameter(blocks, "DYNAMICS", "iseed", args.seed)
    if args.temperature is not None:
        set_parameter(blocks, "ENSEMBLE", "temperature", input_generator.number(args.temperature))
    if args.tau_t is not None and has_thermostat:
        set_parameter(blocks, "ENSEMBLE", "tau_t", input_generator.number(args.tau_t))
    if args.pressure is not None and has_barostat:
        set_parameter(blocks, "ENSEMBLE", "pressure", input_generator.number(args.pressure))
    if args.tau_p is not None and has_barostat:
        set_parameter(blocks, "ENSEMBLE", "tau_p", input_generator.number(args.tau_p))

    cutoff = float_parameter(blocks, "ENERGY", "cutoffdist")
    pairlist = float_parameter(blocks, "ENERGY", "pairlistdist")
    current_skin = pairlist - cutoff
    if current_skin <= 0.0:
        raise RuntimeError("[ENERGY] pairlistdist must be greater than cutoffdist")
    selected_cutoff = args.cutoff if args.cutoff is not None else cutoff
    selected_skin = args.pair_list_skin if args.pair_list_skin is not None else current_skin
    if selected_cutoff <= 0.0:
        raise RuntimeError("[ENERGY] cutoffdist must be positive")
    if variant.family == "CHARMM" and selected_cutoff <= 2.0:
        raise RuntimeError("cutoff is smaller than the CHARMM switching width")
    if args.cutoff is not None or args.pair_list_skin is not None:
        set_parameter(blocks, "ENERGY", "cutoffdist", input_generator.number(selected_cutoff))
        set_parameter(
            blocks, "ENERGY", "pairlistdist",
            input_generator.number(selected_cutoff + selected_skin),
        )
    if has_parameter(blocks, "ENERGY", "switchdist"):
        expected_switch = (selected_cutoff if variant.family == "AMBER"
                           else selected_cutoff - 2.0)
        if args.cutoff is not None:
            set_parameter(
                blocks, "ENERGY", "switchdist",
                input_generator.number(expected_switch),
            )
        elif abs(float_parameter(blocks, "ENERGY", "switchdist") - expected_switch) > 1.0e-9:
            raise RuntimeError("[ENERGY] switchdist does not match native switching policy")

    nsteps = integer_parameter(blocks, "DYNAMICS", "nsteps")
    output_period = integer_parameter(blocks, "DYNAMICS", "eneout_period")
    if nsteps < 1 or output_period < 1:
        raise RuntimeError("nsteps and eneout_period must be positive")
    if integer_parameter(blocks, "DYNAMICS", "iseed") < 1:
        raise RuntimeError("iseed must be positive")
    if nsteps % output_period:
        raise RuntimeError("nsteps must be divisible by eneout_period")
    if has_thermostat:
        thermostat_period = integer_parameter(blocks, "DYNAMICS", "thermostat_period")
        if thermostat_period < 1:
            raise RuntimeError("thermostat_period must be positive")
        if nsteps % thermostat_period:
            raise RuntimeError("NVT/NPT nsteps must be divisible by thermostat_period")
    if has_barostat:
        barostat_period = integer_parameter(blocks, "DYNAMICS", "barostat_period")
        baroscale_period = integer_parameter(blocks, "DYNAMICS", "baroscale_period")
        if barostat_period < 1 or baroscale_period < 1:
            raise RuntimeError("barostat periods must be positive")
        if nsteps % barostat_period:
            raise RuntimeError("NPT nsteps must be divisible by barostat_period")
        if barostat_period % baroscale_period:
            raise RuntimeError("NPT barostat_period must be divisible by baroscale_period")
    if float_parameter(blocks, "ENSEMBLE", "temperature") <= 0.0:
        raise RuntimeError("temperature must be positive")
    if has_thermostat and float_parameter(blocks, "ENSEMBLE", "tau_t") <= 0.0:
        raise RuntimeError("tau_t must be positive")
    if has_barostat:
        if float_parameter(blocks, "ENSEMBLE", "pressure") <= 0.0:
            raise RuntimeError("pressure must be positive")
        if float_parameter(blocks, "ENSEMBLE", "tau_p") <= 0.0:
            raise RuntimeError("tau_p must be positive")


def parse_performance(text):
    """Return the first GENESIS ns/day performance value in text."""
    performance_match = PERFORMANCE_RE.search(text)
    return float(performance_match.group(1)) if performance_match else None


def missing_gpu_route_evidence(text):
    """Return required GENESIS GPU-route labels absent from runtime output."""
    return tuple(label for label, pattern in GPU_ROUTE_PATTERNS
                 if pattern.search(text or "") is None)


def genesis_completion_error(text, expected_steps):
    """Return a scientific-completion error, or ``None`` for a clean run."""
    complete_text = text or ""
    dynamics = complete_text[complete_text.find("[STEP5]"):] \
        if "[STEP5]" in complete_text else complete_text
    if re.search(
            r"(?i)(?:^|\s)(?:nan|[-+]?inf(?:inity)?|\*{3,})(?:\s|$)",
            dynamics,
    ):
        return "non-finite or overflow marker in dynamics output"
    steps = [int(value) for value in re.findall(
        r"^INFO:\s+(\d+)\s+", dynamics, re.MULTILINE,
    )]
    if not steps or steps[-1] != expected_steps:
        return "final dynamics step is %s, expected %d" % (
            steps[-1] if steps else "missing", expected_steps,
        )
    if re.search(r"\[STEP6\]\s*Deallocate Arrays", complete_text) is None:
        return "normal STEP6 deallocation marker is missing"
    if PERFORMANCE_RE.search(complete_text) is None:
        return "performance footer is missing"
    return None


def parse_num_atoms(text):
    """Return the first atom-count value in GENESIS output."""
    for atom_pattern in NUM_ATOMS_PATTERNS:
        atom_match = atom_pattern.search(text)
        if atom_match:
            return int(atom_match.group(1))
    return None


def normalize_csv_key(text):
    """Convert section and parameter names into stable CSV suffixes."""
    return re.sub(r"[^0-9a-zA-Z_]+", "_", text.strip().lower()).strip("_")


def input_options(text):
    """Return section-qualified input parameters from a GENESIS control file."""
    options = {}
    for section_name, section_lines in split_sections(text):
        if not section_name:
            continue
        section_key = normalize_csv_key(section_name)
        line_index = 1
        while line_index < len(section_lines):
            raw_line = section_lines[line_index]
            stripped_line = raw_line.strip()
            if not stripped_line or stripped_line.startswith("#"):
                line_index += 1
                continue
            key_match = KEYLINE_RE.match(raw_line)
            if not key_match:
                raise RuntimeError(
                    "unparsed line in GENESIS section [%s]: %r" %
                    (section_name, stripped_line)
                )

            parameter_key = normalize_csv_key(key_match.group(1))
            value_parts = []
            raw_value = raw_line.split("=", 1)[1]
            while True:
                value_text = raw_value.split("#", 1)[0].strip()
                value_continues = value_text.endswith("\\")
                if value_continues:
                    value_text = value_text[:-1].rstrip()
                if value_text:
                    value_parts.append(value_text)
                line_index += 1
                if not value_continues:
                    break
                if line_index >= len(section_lines):
                    raise RuntimeError(
                        "unterminated continuation for [%s] %s" %
                        (section_name, key_match.group(1))
                    )
                raw_value = section_lines[line_index].strip()
            qualified_key = "input_%s_%s" % (section_key, parameter_key)
            if qualified_key in options:
                raise RuntimeError("duplicate GENESIS option: %s" % qualified_key)
            options[qualified_key] = " ".join(value_parts)
    return options


def canonical_genesis_text(variant, ensemble, dt_fs):
    """Render the authoritative GENESIS base input for one runner cell."""
    period = max(1, int(round(20.0 / dt_fs)))
    protocol = input_generator.Protocol(
        ensemble=ensemble,
        dt_fs=dt_fs,
        nsteps=input_generator.DEFAULT_NSTEPS,
        output_period=input_generator.DEFAULT_OUTPUT_PERIOD,
        thermostat_period=period,
        barostat_period=period,
        baroscale_period=period,
        cutoff_angstrom=(9.0 if variant.family == "AMBER" else 12.0),
        pair_list_skin_angstrom=input_generator.DEFAULT_PAIR_LIST_SKIN,
        temperature_kelvin=input_generator.DEFAULT_TEMPERATURE_K,
        pressure_atm=input_generator.DEFAULT_PRESSURE_ATM,
        random_seed=input_generator.DEFAULT_RANDOM_SEED,
        tau_t_ps=input_generator.DEFAULT_TAU_T_PS,
        tau_p_ps=input_generator.DEFAULT_TAU_P_PS,
        namd_piston_period_fs=input_generator.DEFAULT_NAMD_PISTON_PERIOD_FS,
        namd_piston_decay_fs=input_generator.DEFAULT_NAMD_PISTON_DECAY_FS,
    )
    return input_generator.genesis_input(
        variant, input_generator.SYSTEMS[variant.system],
        variant.asset_for("GENESIS"), protocol,
    )


def validate_canonical_immutable_controls(base_text, variant, ensemble, dt_fs):
    """Reject any base-input drift outside the explicit runner override set."""
    canonical_text = canonical_genesis_text(variant, ensemble, dt_fs)
    actual_sections = tuple(
        normalize_csv_key(name) for name, _lines in split_sections(base_text) if name
    )
    expected_sections = tuple(
        normalize_csv_key(name) for name, _lines in split_sections(canonical_text) if name
    )
    if actual_sections != expected_sections:
        raise RuntimeError(
            "GENESIS section sequence differs (got %r, canonical %r)" %
            (actual_sections, expected_sections)
        )
    actual = input_options(base_text)
    expected = input_options(canonical_text)
    allowed = MUTABLE_BASE_OPTIONS & set(expected)
    actual_immutable = {key: value for key, value in actual.items() if key not in allowed}
    expected_immutable = {key: value for key, value in expected.items() if key not in allowed}
    if actual_immutable == expected_immutable:
        return
    missing = sorted(set(expected_immutable) - set(actual_immutable))
    extra = sorted(set(actual_immutable) - set(expected_immutable))
    changed = sorted(
        key for key in set(actual_immutable) & set(expected_immutable)
        if actual_immutable[key] != expected_immutable[key]
    )
    details = []
    if missing:
        details.append("missing " + ",".join(missing))
    if extra:
        details.append("unexpected " + ",".join(extra))
    if changed:
        key = changed[0]
        details.append("%s=%r, canonical=%r" %
                       (key, actual_immutable[key], expected_immutable[key]))
    raise RuntimeError("immutable GENESIS controls differ (%s)" % "; ".join(details))


def acquire_lock(lock_path, label, poll_seconds=2.0):
    """Acquire the advisory benchmark lock and record its owner."""
    lock_dir = os.path.dirname(os.path.abspath(lock_path))
    if lock_dir:
        os.makedirs(lock_dir, exist_ok=True)
    lock_fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o666)
    while True:
        try:
            fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            os.ftruncate(lock_fd, 0)
            os.lseek(lock_fd, 0, os.SEEK_SET)
            os.write(lock_fd, ("%d %s\n" % (os.getpid(), label)).encode())
            return lock_fd
        except BlockingIOError:
            time.sleep(poll_seconds)


def release_lock(lock_fd):
    """Release and close the advisory benchmark lock file descriptor."""
    try:
        os.ftruncate(lock_fd, 0)
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
    finally:
        os.close(lock_fd)


def safe_tar_member(data_root, system_name, member):
    """Return True when a tar member is safe to extract for a system."""
    member_name = member.name
    expected_prefix = system_name + "/"
    if member_name == system_name:
        if not member.isdir():
            return False
    elif not member_name.startswith(expected_prefix):
        return False

    target_path = os.path.abspath(os.path.join(data_root, member_name))
    system_root = os.path.abspath(os.path.join(data_root, system_name))
    if target_path != system_root and not target_path.startswith(system_root + os.sep):
        return False
    return member.isfile() or member.isdir()


def sha256_file(path):
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as file_obj:
        for chunk in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_command_output(command):
    """Return one-line read-only command output, or ``None`` on failure."""
    try:
        completed = subprocess.run(
            command, cwd=HERE, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL, timeout=5, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    return " | ".join(line.strip() for line in completed.stdout.splitlines() if line.strip())


def runtime_provenance(spdyn_path):
    """Return stable hardware/source facts required to interpret measurements."""
    commit = read_command_output(("git", "rev-parse", "HEAD"))
    tracked_status = read_command_output(("git", "status", "--short", "--untracked-files=no"))
    gpu = read_command_output((
        "nvidia-smi", "--query-gpu=name,driver_version,pci.bus_id",
        "--format=csv,noheader",
    ))
    return {
        "hostname": platform.node() or "unavailable",
        "platform": platform.platform(),
        "processor": platform.processor() or "unavailable",
        "python": platform.python_version(),
        "git_commit": commit or "unavailable",
        "tracked_worktree": ("unavailable" if tracked_status is None
                             else "clean" if not tracked_status else tracked_status),
        "spdyn_sha256": sha256_file(spdyn_path),
        "gpu": gpu or "unavailable",
    }


def system_data_is_current(system_dir, member_hashes, archive_digest):
    """Return True when extracted system data matches every pinned member byte."""
    sentinel_path = os.path.join(system_dir, ARCHIVE_SENTINEL)
    try:
        if open(sentinel_path).read().strip() != archive_digest:
            return False
    except OSError:
        return False

    for member_name, expected_digest in member_hashes.items():
        member_path = os.path.join(DATA, *member_name.split("/"))
        if not os.path.isfile(member_path):
            return False
        if sha256_file(member_path) != expected_digest:
            return False
    return True


def extract_system_archive(archive_path, system_name, members, archive_digest, member_hashes):
    """Extract a validated system archive into data/<system>."""
    temp_root = tempfile.mkdtemp(prefix=".extract-%s-" % system_name, dir=DATA)
    temp_system_dir = os.path.join(temp_root, system_name)
    system_dir = os.path.join(DATA, system_name)
    try:
        with tarfile.open(archive_path, "r:gz") as archive_file:
            archive_file.extractall(temp_root, members)
        if not os.path.isdir(temp_system_dir):
            raise RuntimeError("%s did not contain %s/" % (archive_path, system_name))
        for member_name, expected_digest in member_hashes.items():
            relative = member_name.split("/")[1:]
            extracted = os.path.join(temp_system_dir, *relative)
            if not os.path.isfile(extracted) or sha256_file(extracted) != expected_digest:
                raise RuntimeError("extracted member failed SHA-256 validation: %s" % member_name)
        with open(os.path.join(temp_system_dir, ARCHIVE_SENTINEL), "w") as sentinel_file:
            sentinel_file.write(archive_digest + "\n")
        if os.path.exists(system_dir):
            shutil.rmtree(system_dir)
        os.rename(temp_system_dir, system_dir)
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def ensure_system_data(system_names):
    """Extract data archives for selected systems when needed."""
    extracted_systems = []
    os.makedirs(DATA, exist_ok=True)
    for system_name in system_names:
        system_dir = os.path.join(DATA, system_name)
        archive_path = os.path.join(DATA, system_name + ".tgz")
        if not os.path.isfile(archive_path):
            raise RuntimeError("missing %s and %s" % (system_dir, archive_path))

        pin = input_generator.ARCHIVE_PINS.get(system_name)
        if pin is None:
            raise RuntimeError("no pinned archive provenance for %s" % system_name)
        archive_digest = sha256_file(archive_path)
        if archive_digest != pin["archive_sha256"]:
            raise RuntimeError(
                "%s SHA-256 differs from the pinned native archive" % archive_path
            )
        with tarfile.open(archive_path, "r:gz") as archive_file:
            members = archive_file.getmembers()
            unsafe_members = [
                member.name for member in members
                if not safe_tar_member(DATA, system_name, member)
            ]
            if unsafe_members:
                raise RuntimeError(
                    "unsafe member(s) in %s: %s" %
                    (archive_path, ", ".join(unsafe_members[:5]))
                )
            regular = {member.name: member for member in members if member.isfile()}
            if set(regular) != set(pin["members"]):
                raise RuntimeError(
                    "%s regular member set differs from the pinned registry" % archive_path
                )
            member_hashes = {}
            for member_name in sorted(regular):
                stream = archive_file.extractfile(regular[member_name])
                if stream is None:
                    raise RuntimeError("cannot read %s::%s" % (archive_path, member_name))
                digest = hashlib.sha256()
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
                actual = digest.hexdigest()
                if actual != pin["members"][member_name]:
                    raise RuntimeError("%s::%s failed pinned SHA-256" %
                                       (archive_path, member_name))
                member_hashes[member_name] = actual

        if (os.path.isdir(system_dir)
                and system_data_is_current(system_dir, member_hashes, archive_digest)):
            continue
        extract_system_archive(
            archive_path, system_name, members, archive_digest, member_hashes,
        )
        extracted_systems.append(system_name)
    return extracted_systems


def terminate_process_group(process):
    """Terminate a subprocess process group, escalating to SIGKILL if needed."""
    if process is None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5.0)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def write_spdyn_log(log_path, phase, tag, run_id, input_path, command, environment, return_code,
                    wall_seconds, timeout_seconds, stdout_text, stderr_text):
    """Write one raw GENESIS subprocess log file."""
    if log_path is None:
        return
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as log_file:
        log_file.write("# GENESIS benchmark run log\n")
        log_file.write("# phase: %s\n" % (phase or "unknown"))
        log_file.write("# tag: %s\n" % (tag or "unknown"))
        log_file.write("# run_id: %s\n" % (run_id if run_id is not None else ""))
        log_file.write("# input: %s\n" % input_path)
        log_file.write("# cwd: %s\n" % HERE)
        log_file.write("# command: %s\n" % " ".join(shlex.quote(str(part)) for part in command))
        log_file.write(
            "# env: GENESIS_GPU_PROFILE=%s OMP_NUM_THREADS=%s HWLOC_COMPONENTS=%s\n" %
            (
                environment.get("GENESIS_GPU_PROFILE", ""),
                environment.get("OMP_NUM_THREADS", ""),
                environment.get("HWLOC_COMPONENTS", ""),
            )
        )
        log_file.write("# return_code: %s\n" % return_code)
        log_file.write("# wall_seconds: %.6f\n" % wall_seconds)
        log_file.write("# timeout_seconds: %s\n" % timeout_seconds)
        log_file.write("\n# stdout\n")
        log_file.write(stdout_text or "")
        if stdout_text and not stdout_text.endswith("\n"):
            log_file.write("\n")
        log_file.write("\n# stderr\n")
        log_file.write(stderr_text or "")
        if stderr_text and not stderr_text.endswith("\n"):
            log_file.write("\n")


def run_spdyn(input_path, timeout_seconds, mpi_procs, omp_threads, log_path=None,
              phase=None, tag=None, run_id=None):
    """Run one spdyn command and optionally write its raw log."""
    environment = dict(
        os.environ,
        GENESIS_GPU_PROFILE="0",
        OMP_NUM_THREADS=str(omp_threads),
        HWLOC_COMPONENTS="x86",
    )
    command = ["mpirun", "-np", str(mpi_procs), SPDYN, input_path]
    process = None
    stdout_text = ""
    stderr_text = ""
    return_code = None
    start_time = time.time()
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=environment,
            cwd=HERE,
            start_new_session=True,
        )
        try:
            stdout_text, stderr_text = process.communicate(timeout=timeout_seconds)
            return_code = process.returncode
        except subprocess.TimeoutExpired:
            terminate_process_group(process)
            stdout_text, stderr_text = process.communicate()
            return_code = 124
        wall_seconds = time.time() - start_time
    except KeyboardInterrupt:
        terminate_process_group(process)
        if process is not None:
            try:
                stdout_text, stderr_text = process.communicate(timeout=1.0)
            except subprocess.TimeoutExpired:
                terminate_process_group(process)
                stdout_text, stderr_text = process.communicate()
            except BaseException:
                stdout_text = stdout_text or ""
                stderr_text = stderr_text or ""
        return_code = 130
        wall_seconds = time.time() - start_time
        stderr_text = (stderr_text or "") + "\n[run_benchmark] interrupted by KeyboardInterrupt\n"
        write_spdyn_log(
            log_path, phase, tag, run_id, input_path, command, environment,
            return_code, wall_seconds, timeout_seconds, stdout_text, stderr_text,
        )
        raise

    write_spdyn_log(
        log_path, phase, tag, run_id, input_path, command, environment,
        return_code, wall_seconds, timeout_seconds, stdout_text, stderr_text,
    )
    return stdout_text, stderr_text, return_code, wall_seconds


def build_run_input(base_text, variant, args):
    """Return a validated measurement input without adding GPU parameters."""
    blocks = split_sections(base_text)
    apply_protocol_overrides(blocks, variant, args)
    return join_sections(blocks)


def open_benchmark_log(path):
    """Open the full benchmark log and flush buffered startup messages."""
    global BENCHMARK_LOG_FILE, BENCHMARK_LOG_BUFFER
    os.makedirs(os.path.dirname(path), exist_ok=True)
    BENCHMARK_LOG_FILE = open(path, "w")
    for message in BENCHMARK_LOG_BUFFER:
        BENCHMARK_LOG_FILE.write(message + "\n")
    BENCHMARK_LOG_FILE.flush()
    BENCHMARK_LOG_BUFFER = []


def close_benchmark_log():
    """Close the full benchmark log if it is open."""
    global BENCHMARK_LOG_FILE
    if BENCHMARK_LOG_FILE is not None:
        BENCHMARK_LOG_FILE.close()
        BENCHMARK_LOG_FILE = None


def log(message):
    """Write a progress message to stderr and benchmark.log."""
    sys.stderr.write(message + "\n")
    sys.stderr.flush()
    if BENCHMARK_LOG_FILE is not None:
        BENCHMARK_LOG_FILE.write(message + "\n")
        BENCHMARK_LOG_FILE.flush()
    else:
        BENCHMARK_LOG_BUFFER.append(message)


def emit(message=""):
    """Write a final-report line to stdout and benchmark.log."""
    print(message)
    if BENCHMARK_LOG_FILE is not None:
        BENCHMARK_LOG_FILE.write(message + "\n")
        BENCHMARK_LOG_FILE.flush()


def write_summary_log(path, lines):
    """Write the final aggregate table to summary.log."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as summary_file:
        for line in lines:
            summary_file.write(line + "\n")


def summary_table_lines(results, measure_count, num_steps):
    """Return the final aggregate ns/day table lines."""
    lines = [
        "=== ns/day (mean/median +- std, cv%%) : measure=%d, nsteps=%d ===" %
        (measure_count, num_steps),
        "%-10s %-34s %-5s %-6s %12s %12s %10s %7s  %s" %
        ("system", "variant", "ens", "dt", "mean", "median", "+-std", "cv%", "note"),
        "-" * 122,
    ]
    for result in results:
        lines.append(
            "%-10s %-34s %-5s %-6s %12.2f %12.2f %10.2f %6.1f%%  %s" %
            (
                result["system"], result["variant"], result["ensemble"], result["dt"],
                result["mean"], result["median"], result["stddev"],
                result["cv_percent"], result.get("note", ""),
            )
        )
    return lines


def unique_run_paths(path):
    """Return unused CSV and log-directory paths for a benchmark run."""
    root, extension = os.path.splitext(path)
    log_root = root if extension.lower() == ".csv" else path + ".logs"
    if not os.path.exists(path) and not os.path.exists(log_root):
        return path, log_root
    for suffix_number in range(1, 1000):
        candidate_csv = "%s-%d%s" % (root, suffix_number, extension)
        candidate_log_root = "%s-%d" % (log_root, suffix_number)
        if not os.path.exists(candidate_csv) and not os.path.exists(candidate_log_root):
            return candidate_csv, candidate_log_root
    raise RuntimeError("could not find unused output/log paths for %s" % path)


def run_log_path(log_root, subdir, tag, run_name):
    """Return the raw GENESIS log path for one run."""
    return os.path.join(log_root, subdir, "%s_%s.log" % (tag, run_name))


def positive_int(text):
    """Parse an argparse value as a positive integer."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if value < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return value


def nonnegative_int(text):
    """Parse an argparse value as a nonnegative integer."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if value < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return value


def planned_cells(variants, ensembles, time_steps):
    """Return the complete requested benchmark matrix."""
    return [
        (variant, ensemble, dt_fs)
        for variant in variants
        for ensemble in ensembles
        for dt_fs in time_steps
    ]


def cell_label(cell):
    """Return the input-file stem for a benchmark cell tuple."""
    variant, ensemble, dt_fs = cell
    return input_generator.cell_stem(variant.system, variant.variant, ensemble, dt_fs)


def run_cell(variant, ensemble, dt_fs, args):
    """Run warmup and measurement repetitions for one benchmark cell."""
    system_name = variant.system
    tag = cell_label((variant, ensemble, dt_fs))
    base_path = os.path.join(args.input_root, tag + ".inp")
    if not os.path.isfile(base_path):
        raise RuntimeError("missing preflighted base input %s" % base_path)

    base_text = open(base_path).read()
    note = ""
    num_atoms = None

    run_input = os.path.join(args.input_log_dir, tag + ".run.inp")
    run_text = build_run_input(base_text, variant, args)
    open(run_input, "w").write(run_text)
    actual_baroscale_period = get_parameter_value(run_text, "DYNAMICS", "baroscale_period")
    options = input_options(run_text)

    for warmup_run in range(args.warmup):
        log("  [warmup %d/%d] %s" % (warmup_run + 1, args.warmup, tag))
        warmup_log = run_log_path(args.log_root, "production", tag, "warmup_%d" % (warmup_run + 1))
        stdout_text, stderr_text, return_code, wall_seconds = run_spdyn(
            run_input, args.timeout, args.mpi_procs, args.omp_threads,
            warmup_log, "production-warmup", tag, warmup_run + 1,
        )
        del wall_seconds
        num_atoms = parse_num_atoms((stdout_text or "") + "\n" + (stderr_text or "")) or num_atoms
        if return_code != 0:
            log("  [warmup FAILED rc=%d] %s" % (return_code, (stderr_text or stdout_text)[-400:]))
            return None
        completion_error = genesis_completion_error(
            (stdout_text or "") + "\n" + (stderr_text or ""), args.num_steps,
        )
        if completion_error:
            log("  [warmup FAILED: %s]" % completion_error)
            return None
        missing_gpu = missing_gpu_route_evidence(
            (stdout_text or "") + "\n" + (stderr_text or ""),
        )
        if missing_gpu:
            log("  [warmup FAILED: missing GPU route evidence: %s]" %
                ", ".join(missing_gpu))
            return None
        if num_atoms is not None and num_atoms != variant.atom_count:
            log("  [warmup FAILED atom count %d != manifest %d]" %
                (num_atoms, variant.atom_count))
            return None

    runs = []
    for measure_run in range(args.measure):
        production_log = run_log_path(args.log_root, "production", tag, measure_run + 1)
        stdout_text, stderr_text, return_code, wall_seconds = run_spdyn(
            run_input, args.timeout, args.mpi_procs, args.omp_threads,
            production_log, "production-measure", tag, measure_run + 1,
        )
        num_atoms = parse_num_atoms((stdout_text or "") + "\n" + (stderr_text or "")) or num_atoms
        if return_code != 0:
            log("  [measure FAILED rc=%d] %s" % (return_code, (stderr_text or stdout_text)[-400:]))
            return None
        completion_error = genesis_completion_error(
            (stdout_text or "") + "\n" + (stderr_text or ""), args.num_steps,
        )
        if completion_error:
            log("  [measure FAILED: %s]" % completion_error)
            return None
        missing_gpu = missing_gpu_route_evidence(
            (stdout_text or "") + "\n" + (stderr_text or ""),
        )
        if missing_gpu:
            log("  [measure FAILED: missing GPU route evidence: %s]" %
                ", ".join(missing_gpu))
            return None
        if num_atoms is not None and num_atoms != variant.atom_count:
            log("  [measure FAILED atom count %d != manifest %d]" %
                (num_atoms, variant.atom_count))
            return None
        ns_per_day = parse_performance(stdout_text)
        if ns_per_day is None:
            log("  [measure FAILED: no [PERFORMANCE] parsed] %s" % ((stderr_text or stdout_text)[-400:]))
            return None
        runs.append(dict(
            run_id=measure_run + 1,
            ns_per_day=ns_per_day,
            wall_seconds=wall_seconds,
            production_log=os.path.relpath(production_log, HERE),
        ))
        log("  [measure %d/%d] %s  ns/day=%.2f  (wall %.1fs)" %
            (measure_run + 1, args.measure, tag, ns_per_day, wall_seconds))

    if num_atoms is None:
        log("  [measure FAILED: no atom count parsed]")
        return None
    performance_values = [run["ns_per_day"] for run in runs]
    mean_value = statistics.mean(performance_values)
    median_value = statistics.median(performance_values)
    std_value = statistics.stdev(performance_values) if len(performance_values) > 1 else 0.0
    cv_value = (std_value / mean_value * 100.0) if mean_value else 0.0
    return dict(
        system=system_name,
        variant=variant.variant,
        forcefield=variant.forcefield,
        water_model=variant.water_model,
        solvent=variant.solvent,
        ensemble=ensemble,
        dt=input_generator.dt_label(dt_fs),
        mean=mean_value,
        median=median_value,
        stddev=std_value,
        cv_percent=cv_value,
        measure_count=len(performance_values),
        runs=runs,
        note=note,
        num_atoms=num_atoms,
        input_options=options,
        baroscale_period=actual_baroscale_period,
        mpi_procs=args.mpi_procs,
        omp_threads=args.omp_threads,
    )


def csv_fieldnames(input_columns):
    """Return the CSV header for benchmark result rows."""
    return [
        "system", "variant", "forcefield", "water_model", "solvent",
        "ensemble", "dt", "run_id", "ns_per_day", "wall_seconds",
        "ns_per_day_mean", "ns_per_day_median", "ns_per_day_std", "cv_pct",
        "n_measure", "mpi_procs", "omp_threads", "num_atoms",
        "baroscale_period", "note", "production_log", "input_options_json",
    ] + list(input_columns)


def planned_input_columns(cells, args):
    """Return input-option CSV columns for the selected benchmark matrix."""
    columns = set()
    for cell in cells:
        tag = cell_label(cell)
        base_path = os.path.join(args.input_root, tag + ".inp")
        if not os.path.isfile(base_path):
            raise RuntimeError("missing preflighted base input %s" % base_path)
        base_text = open(base_path).read()
        run_text = build_run_input(base_text, cell[0], args)
        columns.update(input_options(run_text))
    return sorted(columns)


def result_csv_rows(result, input_columns):
    """Return CSV row dictionaries for one completed benchmark cell."""
    options = result["input_options"]
    missing_columns = sorted(set(options) - set(input_columns))
    if missing_columns:
        raise RuntimeError("CSV header is missing input option column(s): %s" %
                           ",".join(missing_columns))

    base_row = {
        "system": result["system"],
        "variant": result["variant"],
        "forcefield": result["forcefield"],
        "water_model": result["water_model"],
        "solvent": result["solvent"],
        "ensemble": result["ensemble"],
        "dt": result["dt"],
        "ns_per_day_mean": "%.3f" % result["mean"],
        "ns_per_day_median": "%.3f" % result["median"],
        "ns_per_day_std": "%.3f" % result["stddev"],
        "cv_pct": "%.2f" % result["cv_percent"],
        "n_measure": result["measure_count"],
        "mpi_procs": result["mpi_procs"],
        "omp_threads": result["omp_threads"],
        "num_atoms": result["num_atoms"] or "",
        "baroscale_period": result["baroscale_period"] or "",
        "note": result.get("note", ""),
        "input_options_json": json.dumps(options, sort_keys=True, separators=(",", ":")),
    }
    base_row.update(options)

    rows = []
    for run in result["runs"]:
        row = dict(base_row)
        row.update({
            "run_id": run["run_id"],
            "ns_per_day": "%.3f" % run["ns_per_day"],
            "wall_seconds": "%.3f" % run["wall_seconds"],
            "production_log": run["production_log"],
        })
        rows.append(row)
    return rows


def create_results_csv(csv_path, input_columns):
    """Create a result CSV with a header row."""
    with open(csv_path, "x", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_fieldnames(input_columns))
        writer.writeheader()


def append_result_csv(csv_path, result, input_columns):
    """Append one completed cell to the result CSV."""
    with open(csv_path, "a", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_fieldnames(input_columns), extrasaction="raise")
        writer.writerows(result_csv_rows(result, input_columns))
        csv_file.flush()
        os.fsync(csv_file.fileno())


def build_parser():
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--systems", help="comma list; default all systems with available native models")
    parser.add_argument("--ensembles", help="comma list of nve,nvt,npt (default all valid)")
    parser.add_argument("--dt", default="2,4", help="comma list of time steps in fs")
    parser.add_argument("--warmup", type=nonnegative_int, default=2,
                        help="discarded warm-up runs (default 2)")
    parser.add_argument("--measure", type=positive_int, default=10,
                        help="timed runs (default 10)")
    parser.add_argument("--nsteps", dest="num_steps", metavar="NSTEPS", type=positive_int, default=100000,
                        help="measurement-run nsteps (default 100000)")
    parser.add_argument("--eneout-period", "--output-period", dest="eneout_period",
                        type=positive_int, default=1000,
                        help="DYNAMICS eneout_period for measurement inputs (default 1000)")
    parser.add_argument("--thermostat-period", type=positive_int,
                        help="DYNAMICS thermostat_period override")
    parser.add_argument("--barostat-period", type=positive_int,
                        help="DYNAMICS barostat_period override")
    parser.add_argument("--mpi-procs", type=positive_int, default=1,
                        help="MPI process count for each spdyn launch (default 1)")
    parser.add_argument("--omp-threads", type=positive_int, default=1,
                        help="OMP_NUM_THREADS for each spdyn launch (default 1)")
    parser.add_argument("--baroscale-period", type=positive_int, default=None,
                        help="DYNAMICS baroscale_period override for NPT work inputs")
    parser.add_argument("--tau-t", type=input_generator.finite_positive,
                        help="ENSEMBLE tau_t override in ps")
    parser.add_argument("--tau-p", type=input_generator.finite_positive,
                        help="ENSEMBLE tau_p override in ps")
    parser.add_argument("--temperature", type=input_generator.finite_positive,
                        help="ENSEMBLE temperature override in kelvin")
    parser.add_argument("--pressure", type=input_generator.finite_positive,
                        help="ENSEMBLE pressure override in atm")
    parser.add_argument("--seed", type=positive_int,
                        help="DYNAMICS iseed override")
    parser.add_argument("--cutoff", type=input_generator.finite_positive,
                        help="ENERGY cutoffdist override in Angstrom")
    parser.add_argument("--pair-list-skin", "--pairlist-skin", dest="pair_list_skin",
                        type=input_generator.finite_positive,
                        help="pair-list skin override in Angstrom")
    parser.add_argument("--input-root", default=INPUTS,
                        help="directory containing generated GENESIS inputs")
    parser.add_argument("--timeout", type=positive_int, default=7200,
                        help="per-run timeout seconds (default 7200)")
    parser.add_argument("--lock", default="/tmp/bench.lock",
                        help="advisory serialisation lockfile")
    parser.add_argument("--allow-failures", action="store_true",
                        help="write a partial CSV and exit 0 even if selected cells fail")
    parser.add_argument("--timestamp", default=None,
                        help="results filename stamp (else auto)")
    parser.add_argument("--out", default=None,
                        help="explicit CSV output path (mutually exclusive with --timestamp)")
    parser.add_argument("--spdyn", default=SPDYN,
                        help="path to the spdyn binary")
    return parser


def parse_time_steps(time_steps_text, parser):
    """Return unique, finite selected time steps in femtoseconds."""
    time_steps = []
    for time_step_text in time_steps_text.split(","):
        normalized = time_step_text.strip().lower()
        if normalized.endswith("fs"):
            normalized = normalized[:-2]
        if not normalized:
            continue
        try:
            value = float(normalized)
        except ValueError:
            parser.error("--dt contains a nonnumeric value %r" % time_step_text)
        if not math.isfinite(value) or value <= 0.0 or value > 4.0:
            parser.error("--dt values must be finite and in (0, 4] fs")
        if value not in time_steps:
            time_steps.append(value)
    if not time_steps:
        parser.error("--dt must contain at least one value")
    labels = {}
    for value in time_steps:
        label = input_generator.dt_label(value)
        if label in labels and value != labels[label]:
            parser.error("--dt values %r and %r collide as %s" %
                         (labels[label], value, label))
        labels[label] = value
    return time_steps


def validate_unique_cell_paths(cells, input_root, parser):
    """Reject any matrix aliases that map to the same canonical input path."""
    planned = {}
    for cell in cells:
        path = os.path.realpath(os.path.join(input_root, cell_label(cell) + ".inp"))
        if path in planned and planned[path] != cell:
            parser.error("benchmark cells %r and %r collide at %s" %
                         (planned[path], cell, path))
        planned[path] = cell


def canonical_requested_csv(args, input_root, parser):
    """Resolve and validate the requested result path before acquiring the lock."""
    if args.out is not None and args.timestamp is not None:
        parser.error("--out and --timestamp are mutually exclusive")
    if args.out is not None:
        path = os.path.realpath(os.path.expanduser(args.out))
    else:
        stamp = args.timestamp or datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", stamp)
                or stamp in (".", "..")):
            parser.error("--timestamp must be a filename-safe label")
        path = os.path.realpath(os.path.join(RESULTS, stamp + ".csv"))
    input_root = os.path.realpath(input_root)
    try:
        inside_inputs = os.path.commonpath((path, input_root)) == input_root
    except ValueError:
        inside_inputs = False
    if inside_inputs:
        parser.error("result output must not be inside the generated input tree")
    return path


def parse_ensembles(ensembles_text, parser):
    """Return unique, validated ensemble names."""
    ensembles = input_generator.comma_values(
        ensembles_text or ",".join(ALL_ENSEMBLES), lower=True,
    )
    input_generator.validate_choices(ensembles, ALL_ENSEMBLES, "--ensembles", parser)
    if not ensembles:
        parser.error("--ensembles must contain at least one value")
    return ensembles


def select_variants(args, parser):
    """Resolve the same native registry used by input generation."""
    try:
        loaded, missing = input_generator.available_variants()
        variants, strict = input_generator.select_variants(args, parser, loaded, missing)
    except input_generator.GenerationError as error:
        parser.error(str(error))
    del strict
    selected = []
    for variant in variants:
        genesis_asset = variant.asset_for("GENESIS")
        if genesis_asset is None:
            if args.systems:
                parser.error("%s/%s has no GENESIS asset" % (variant.system, variant.variant))
            continue
        selected.append(variant)
    if not selected:
        parser.error("no selected variant has a GENESIS asset")
    return tuple(selected)


def preflight_inputs(cells, args, parser):
    """Validate every selected generated control before acquiring the run lock."""
    missing = []
    errors = []
    for cell in cells:
        variant, ensemble, dt_fs = cell
        path = os.path.join(args.input_root, cell_label(cell) + ".inp")
        if not os.path.isfile(path):
            missing.append(path)
            continue
        try:
            with open(path, encoding="utf-8") as stream:
                base_text = stream.read()
            required_headers = (
                "# %s" % input_generator.GENERATED_MARKER,
                "# benchmark-system: %s" % variant.system,
                "# native-model: %s" % variant.variant,
                "# forcefield: %s" % variant.forcefield,
                "# water-model: %s" % variant.water_model,
                "# solvent: %s" % variant.solvent,
            )
            absent = [header for header in required_headers if header not in base_text]
            if absent:
                raise RuntimeError("missing provenance header(s): %s" % ", ".join(absent))
            validate_canonical_immutable_controls(
                base_text, variant, ensemble, dt_fs,
            )
            input_ensemble = get_parameter_value(base_text, "ENSEMBLE", "ensemble")
            if input_ensemble is None or input_ensemble.lower() != ensemble:
                raise RuntimeError("[ENSEMBLE] ensemble does not match filename")
            timestep = get_parameter_value(base_text, "DYNAMICS", "timestep")
            if timestep is None or abs(float(timestep) * 1000.0 - dt_fs) > 1.0e-9:
                raise RuntimeError("[DYNAMICS] timestep does not match filename")
            build_run_input(base_text, variant, args)
        except (OSError, ValueError, RuntimeError) as error:
            errors.append("%s: %s" % (path, error))
    if missing:
        parser.error(
            "missing %d generated GENESIS input(s), first: %s; run generate_inputs.py "
            "with matching filters" % (len(missing), missing[0])
        )
    if errors:
        parser.error("input preflight failed: %s" % errors[0])


def main():
    """Run the benchmark command-line program."""
    global SPDYN

    parser = build_parser()
    args = parser.parse_args()
    args.input_root = os.path.realpath(os.path.expanduser(args.input_root))
    requested_csv_path = canonical_requested_csv(args, args.input_root, parser)
    args.lock = os.path.realpath(os.path.expanduser(args.lock))
    variants = select_variants(args, parser)
    ensembles = parse_ensembles(args.ensembles, parser)
    time_steps = parse_time_steps(args.dt, parser)
    cells = planned_cells(variants, ensembles, time_steps)
    if not cells:
        parser.error("selection produced no benchmark cells")
    validate_unique_cell_paths(cells, args.input_root, parser)
    preflight_inputs(cells, args, parser)
    system_names = list(dict.fromkeys(variant.system for variant in variants))

    SPDYN = os.path.realpath(os.path.expanduser(args.spdyn))
    if not (os.path.isfile(SPDYN) and os.access(SPDYN, os.X_OK)):
        sys.exit("error: spdyn not found or not executable: %s" % SPDYN)

    results = []
    failures = []
    csv_path = None
    log_root = None
    log("waiting for benchmark lock: %s" % args.lock)
    lock_fd = acquire_lock(args.lock, "run_benchmark")
    try:
        try:
            archive_native_systems = list(dict.fromkeys(
                variant.system for variant in variants
                if variant.asset_for("GENESIS").archive_native
            ))
            extracted_systems = (ensure_system_data(archive_native_systems)
                                 if archive_native_systems else [])
            csv_path, log_root = unique_run_paths(requested_csv_path)
            csv_path = os.path.realpath(csv_path)
            log_root = os.path.realpath(log_root)
            if csv_path == log_root:
                raise RuntimeError("CSV and log outputs resolve to the same path")
            if os.path.commonpath((csv_path, log_root)) == log_root:
                raise RuntimeError("CSV output must not be inside its log directory")
            if args.lock == csv_path or args.lock == log_root:
                raise RuntimeError("benchmark lock collides with a result output path")
            try:
                log_inputs_overlap = (
                    os.path.commonpath((log_root, args.input_root)) in
                    (log_root, args.input_root)
                )
            except ValueError:
                log_inputs_overlap = False
            if log_inputs_overlap:
                raise RuntimeError("benchmark logs and generated inputs must be disjoint trees")
            args.log_root = log_root
            args.input_log_dir = os.path.join(log_root, "inputs")
            os.makedirs(os.path.join(log_root, "production"), exist_ok=True)
            os.makedirs(args.input_log_dir, exist_ok=True)
            open_benchmark_log(os.path.join(log_root, "benchmark.log"))

            input_columns = planned_input_columns(cells, args)
            create_results_csv(csv_path, input_columns)

            log("=== GENESIS GPU benchmark ===")
            log("spdyn      : %s" % SPDYN)
            provenance = runtime_provenance(SPDYN)
            for key in ("spdyn_sha256", "git_commit", "tracked_worktree", "hostname",
                        "platform", "processor", "python", "gpu"):
                log("%-11s: %s" % (key, provenance[key]))
            log("systems    : %s" % ",".join(system_names))
            log("models     : %s" % ",".join(
                "%s/%s" % (variant.system, variant.variant) for variant in variants
            ))
            log("ensembles  : %s" % ",".join(ensembles))
            log("dt         : %s" % ",".join(input_generator.dt_label(value) for value in time_steps))
            log("warmup/measure: %d/%d   nsteps: %d" %
                (args.warmup, args.measure, args.num_steps))
            log("eneout     : %d" % args.eneout_period)
            log("mpi/omp    : %d/%d" % (args.mpi_procs, args.omp_threads))
            log("baroscale : %s" % (args.baroscale_period if args.baroscale_period is not None else "input"))
            log("lock       : %s" % args.lock)
            log("csv        : %s" % csv_path)
            log("logs       : %s" % log_root)
            log("benchmark  : %s" % os.path.join(log_root, "benchmark.log"))
            log("summary    : %s" % os.path.join(log_root, "summary.log"))
            if csv_path != requested_csv_path:
                log("csv note   : requested output existed; using non-overwriting suffix")
            for system_name in extracted_systems:
                log("data       : extracted data/%s.tgz -> data/%s/" % (system_name, system_name))
            log("")

            for variant, ensemble, dt_fs in cells:
                result = run_cell(variant, ensemble, dt_fs, args)
                if result:
                    results.append(result)
                    append_result_csv(csv_path, result, input_columns)
                else:
                    failures.append(cell_label((variant, ensemble, dt_fs)))
        finally:
            release_lock(lock_fd)

        table_lines = summary_table_lines(results, args.measure, args.num_steps)
        emit("")
        for line in table_lines:
            emit(line)
        write_summary_log(os.path.join(log_root, "summary.log"), table_lines)
        emit("")
        emit("CSV written to %s" % csv_path)
        emit("Logs written to %s" % log_root)
        if failures:
            log("")
            log("FAILED cells: %s" % ",".join(failures))
            if not args.allow_failures:
                sys.exit(1)
    finally:
        close_benchmark_log()


if __name__ == "__main__":
    main()
