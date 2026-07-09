#!/usr/bin/env python3
"""Run GENESIS spdyn GPU benchmarks.

For each selected system, ensemble, and time step, the driver runs an optional
autotuning pass, writes a pinned input, warms up, measures production runs, and
records raw logs plus a row-per-run CSV.
"""

import argparse
import csv
import datetime
import fcntl
import hashlib
import json
import os
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


HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
INPUTS = os.path.join(HERE, "inputs")
RESULTS = os.path.join(HERE, "results")
SPDYN = os.path.normpath(os.path.join(HERE, "..", "src", "spdyn_singlempi", "spdyn"))

ALL_SYSTEMS = ["dhfr", "apoa1", "uun", "factorix", "bpti", "dppc", "ake", "stmv"]
ALL_ENSEMBLES = ["nve", "nvt", "npt"]
ALL_TIME_STEPS = ["2fs", "4fs"]
ARCHIVE_SENTINEL = ".archive_sha256"
KERNEL_KEYS = [
    "kernel_pme_spread_block",
    "kernel_pme_influence_block",
    "kernel_bonded_block",
    "kernel_constraints_block",
    "kernel_nonbond_inter_minblocks",
    "kernel_nonbond_intra_minblocks",
]

PERFORMANCE_RE = re.compile(r"\[PERFORMANCE\]\s*(?:performance:\s*)?([0-9.]+)\s*ns/day")
NUM_ATOMS_PATTERNS = [
    re.compile(r"^\s*num_atoms\s*=\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
    re.compile(r"^\s*natom\s*=\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
    re.compile(r"^\s*number\s+of\s+atoms\s*[:=]\s*([0-9]+)\b", re.IGNORECASE | re.MULTILINE),
]
SECTION_RE = re.compile(r"^\s*\[([A-Za-z0-9_]+)\]\s*$")
KEYLINE_RE = re.compile(r"^\s*([A-Za-z_]\w*)\s*=")

BENCHMARK_LOG_FILE = None
BENCHMARK_LOG_BUFFER = []


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
    return "\n".join(output_lines).rstrip("\n") + "\n"


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


def set_period_overrides(blocks, nbupdate_period, baroscale_period):
    """Apply optional update-period overrides to a control file."""
    if nbupdate_period is not None:
        set_parameter(blocks, "DYNAMICS", "nbupdate_period", nbupdate_period)
    if baroscale_period is not None and has_parameter(blocks, "DYNAMICS", "barostat_period"):
        set_parameter(blocks, "DYNAMICS", "baroscale_period", baroscale_period)


def parse_autotune_report(text):
    """Extract tuned kernel, cell-size, and neighbor-list values."""
    tuned_values = {}
    report_lines = text.splitlines()

    for line_index, raw_line in enumerate(report_lines):
        if "ms/window" not in raw_line or "(selected)" not in raw_line:
            continue

        pairlistdist = None
        nbupdate_period = None
        for candidate_line in report_lines[line_index + 1:]:
            if (
                "ms/window" in candidate_line
                or "Cell size autotune" in candidate_line
                or "Selected" in candidate_line
            ):
                break
            if pairlistdist is None:
                pairlistdist_match = re.search(r"pairlistdist\s*=\s*([0-9.]+)", candidate_line)
                if pairlistdist_match:
                    pairlistdist = float(pairlistdist_match.group(1))
            if nbupdate_period is None:
                nbupdate_match = re.search(r"nbupdate_period\s*=\s*([0-9]+)", candidate_line)
                if nbupdate_match:
                    nbupdate_period = int(nbupdate_match.group(1))
            if pairlistdist is not None and nbupdate_period is not None:
                break
        if pairlistdist is not None:
            tuned_values["pairlistdist"] = pairlistdist
        if nbupdate_period is not None:
            tuned_values["nbupdate_period"] = nbupdate_period
        break

    for line_index, raw_line in enumerate(report_lines):
        if "Selected configuration:" not in raw_line:
            continue
        for cell_line in report_lines[line_index + 1: line_index + 4]:
            cell_size_match = re.search(r"cell_size\s*=\s*([0-9.]+)", cell_line)
            if cell_size_match:
                tuned_values["cell_size"] = float(cell_size_match.group(1))
                break
        break

    kernel_values = {}
    for raw_line in report_lines:
        kernel_match = re.match(r"\s*(kernel_\w+)\s*=\s*([0-9]+)", raw_line)
        if kernel_match and kernel_match.group(1) in KERNEL_KEYS:
            kernel_values[kernel_match.group(1)] = int(kernel_match.group(2))
    if kernel_values:
        tuned_values["kernel"] = kernel_values
    return tuned_values


def parse_performance(text):
    """Return the first GENESIS ns/day performance value in text."""
    performance_match = PERFORMANCE_RE.search(text)
    return float(performance_match.group(1)) if performance_match else None


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
                line_index += 1
                continue

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
                if not value_continues or line_index >= len(section_lines):
                    break
                raw_value = section_lines[line_index].strip()
            options["input_%s_%s" % (section_key, parameter_key)] = " ".join(value_parts)
    return options


def validate_tuned_values(tuned_values, tune_set):
    """Return requested autotune values missing from a parsed report."""
    missing_values = []
    if "kernel" in tune_set:
        kernel_values = tuned_values.get("kernel") or {}
        missing_values.extend(kernel_key for kernel_key in KERNEL_KEYS if kernel_key not in kernel_values)
    if "cell" in tune_set and "cell_size" not in tuned_values:
        missing_values.append("cell_size")
    if "nblist" in tune_set:
        if "pairlistdist" not in tuned_values:
            missing_values.append("pairlistdist")
        if "nbupdate_period" not in tuned_values:
            missing_values.append("nbupdate_period")
    return missing_values


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


def system_data_is_current(system_dir, regular_members, archive_digest):
    """Return True when extracted system data matches its archive."""
    sentinel_path = os.path.join(system_dir, ARCHIVE_SENTINEL)
    try:
        if open(sentinel_path).read().strip() != archive_digest:
            return False
    except OSError:
        return False

    for member in regular_members:
        member_path = os.path.join(DATA, member.name)
        if not os.path.isfile(member_path):
            return False
        if os.path.getsize(member_path) != member.size:
            return False
    return True


def extract_system_archive(archive_path, system_name, members, archive_digest):
    """Extract a validated system archive into data/<system>."""
    temp_root = tempfile.mkdtemp(prefix=".extract-%s-" % system_name, dir=DATA)
    temp_system_dir = os.path.join(temp_root, system_name)
    system_dir = os.path.join(DATA, system_name)
    try:
        with tarfile.open(archive_path, "r:gz") as archive_file:
            archive_file.extractall(temp_root, members)
        if not os.path.isdir(temp_system_dir):
            raise RuntimeError("%s did not contain %s/" % (archive_path, system_name))
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

        archive_digest = sha256_file(archive_path)
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

        regular_members = [member for member in members if member.isfile()]
        if os.path.isdir(system_dir) and system_data_is_current(system_dir, regular_members, archive_digest):
            continue
        extract_system_archive(archive_path, system_name, members, archive_digest)
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


def build_tune_input(base_text, tune_set, num_steps, eneout_period, nbupdate_period, baroscale_period):
    """Return a tuning input for one benchmark cell."""
    blocks = split_sections(base_text)
    set_run_window(blocks, num_steps, eneout_period)
    set_period_overrides(blocks, nbupdate_period, baroscale_period)
    if "cell" in tune_set:
        set_parameter(blocks, "ENERGY", "cell_size_autotune", "YES")
    if "nblist" in tune_set:
        set_parameter(blocks, "DYNAMICS", "nbupdate_autotune", "YES")
    if "kernel" in tune_set:
        set_parameter(blocks, "GPU", "kernel_autotune", "YES")
    return join_sections(blocks)


def build_pinned_input(base_text, tuned_values, tune_set, num_steps, eneout_period,
                       nbupdate_period, baroscale_period):
    """Return a measurement input with selected autotune values pinned."""
    blocks = split_sections(base_text)
    set_run_window(blocks, num_steps, eneout_period)
    set_parameter(blocks, "ENERGY", "cell_size_autotune", "NO")
    set_parameter(blocks, "DYNAMICS", "nbupdate_autotune", "NO")
    set_parameter(blocks, "GPU", "kernel_autotune", "NO")

    if "kernel" in tune_set and tuned_values.get("kernel"):
        for kernel_key in KERNEL_KEYS:
            if kernel_key in tuned_values["kernel"]:
                set_parameter(blocks, "GPU", kernel_key, tuned_values["kernel"][kernel_key])
    if "cell" in tune_set and "cell_size" in tuned_values:
        set_parameter(blocks, "ENERGY", "cell_size", "%.3f" % tuned_values["cell_size"])
    if "nblist" in tune_set:
        if "pairlistdist" in tuned_values:
            set_parameter(blocks, "ENERGY", "pairlistdist", "%.3f" % tuned_values["pairlistdist"])
        if "nbupdate_period" in tuned_values:
            set_parameter(blocks, "DYNAMICS", "nbupdate_period", tuned_values["nbupdate_period"])
    set_period_overrides(blocks, nbupdate_period, baroscale_period)
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


def summary_table_lines(results, tune_set, measure_count, num_steps):
    """Return the final aggregate ns/day table lines."""
    lines = [
        "=== ns/day (mean/median +- std, cv%%) : tuners=%s, measure=%d, nsteps=%d ===" %
        (",".join(sorted(tune_set)) or "none", measure_count, num_steps),
        "%-10s %-5s %-4s %12s %12s %10s %7s  %s" %
        ("system", "ens", "dt", "mean", "median", "+-std", "cv%", "note"),
        "-" * 82,
    ]
    for result in results:
        lines.append("%-10s %-5s %-4s %12.2f %12.2f %10.2f %6.1f%%  %s" %
                     (
                         result["system"], result["ensemble"], result["dt"],
                         result["mean"], result["median"], result["std"],
                         result["cv"], result.get("note", ""),
                     ))
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


def validate_eneout_period(args, parser):
    """Reject run windows that are not multiples of eneout_period."""
    if args.num_steps % args.eneout_period != 0:
        parser.error("--nsteps (%d) must be a multiple of --eneout-period (%d)" %
                     (args.num_steps, args.eneout_period))
    if args.tune_num_steps % args.eneout_period != 0:
        parser.error("--tune-nsteps (%d) must be a multiple of --eneout-period (%d)" %
                     (args.tune_num_steps, args.eneout_period))


def run_cell(system_name, ensemble, time_step, args):
    """Run tuning, warmup, and measurement for one benchmark cell."""
    tag = "%s_%s_%s" % (system_name, ensemble, time_step)
    base_path = os.path.join(INPUTS, tag + ".inp")
    if not os.path.isfile(base_path):
        log("  [skip] no base input %s" % base_path)
        return None

    base_text = open(base_path).read()
    tuned_values = {}
    note = ""
    num_atoms = None
    autotune_log = ""

    if args.tune_set:
        tune_input = os.path.join(args.input_log_dir, tag + ".tune.inp")
        open(tune_input, "w").write(build_tune_input(
            base_text, args.tune_set, args.tune_num_steps, args.eneout_period,
            args.nbupdate_period, args.baroscale_period,
        ))
        tune_log = run_log_path(args.log_root, "autotune", tag, 1)
        autotune_log = os.path.relpath(tune_log, HERE)
        log("  [tune] %s  tuners=%s  nsteps=%d" %
            (tag, ",".join(sorted(args.tune_set)), args.tune_num_steps))
        stdout_text, stderr_text, return_code, wall_seconds = run_spdyn(
            tune_input, args.timeout, args.mpi_procs, args.omp_threads,
            tune_log, "autotune", tag, 1,
        )
        del wall_seconds
        combined_output = (stdout_text or "") + "\n" + (stderr_text or "")
        num_atoms = parse_num_atoms(combined_output) or num_atoms
        if return_code != 0:
            reason = "cell overflow" if "cell too large" in combined_output else ("rc=%d" % return_code)
            if args.tune_set <= {"kernel"} or "cell overflow" not in reason:
                log("  [tune FAILED: %s] %s" % (reason, (stderr_text or stdout_text)[-400:]))
                return None
            log("  [tune FAILED: %s -> measuring with DEFAULTS] %s" % (reason, tag))
            note = "tune-failed(%s)" % reason
            tuned_values = {}
        else:
            tuned_values = parse_autotune_report(stdout_text)
            missing_values = validate_tuned_values(tuned_values, args.tune_set)
            if missing_values:
                log("  [tune FAILED: missing parsed value(s): %s] %s" %
                    (",".join(missing_values), tag))
                return None
            tuned_message = ", ".join(
                "%s=%s" % (key, value) for key, value in tuned_values.items()
                if key != "kernel"
            )
            if "kernel" in tuned_values:
                tuned_message += " kernel=" + str(tuned_values["kernel"])
            log("    tuned: " + (tuned_message or "(nothing parsed)"))

    pinned_input = os.path.join(args.input_log_dir, tag + ".pinned.inp")
    pinned_text = build_pinned_input(
        base_text, tuned_values, args.tune_set, args.num_steps, args.eneout_period,
        args.nbupdate_period, args.baroscale_period,
    )
    open(pinned_input, "w").write(pinned_text)
    actual_nbupdate_period = get_parameter_value(pinned_text, "DYNAMICS", "nbupdate_period")
    actual_baroscale_period = get_parameter_value(pinned_text, "DYNAMICS", "baroscale_period")
    options = input_options(pinned_text)

    for warmup_run in range(args.warmup):
        log("  [warmup %d/%d] %s" % (warmup_run + 1, args.warmup, tag))
        warmup_log = run_log_path(args.log_root, "production", tag, "warmup_%d" % (warmup_run + 1))
        stdout_text, stderr_text, return_code, wall_seconds = run_spdyn(
            pinned_input, args.timeout, args.mpi_procs, args.omp_threads,
            warmup_log, "production-warmup", tag, warmup_run + 1,
        )
        del wall_seconds
        num_atoms = parse_num_atoms((stdout_text or "") + "\n" + (stderr_text or "")) or num_atoms
        if return_code != 0:
            log("  [warmup FAILED rc=%d] %s" % (return_code, (stderr_text or stdout_text)[-400:]))
            return None

    runs = []
    for measure_run in range(args.measure):
        production_log = run_log_path(args.log_root, "production", tag, measure_run + 1)
        stdout_text, stderr_text, return_code, wall_seconds = run_spdyn(
            pinned_input, args.timeout, args.mpi_procs, args.omp_threads,
            production_log, "production-measure", tag, measure_run + 1,
        )
        num_atoms = parse_num_atoms((stdout_text or "") + "\n" + (stderr_text or "")) or num_atoms
        if return_code != 0:
            log("  [measure FAILED rc=%d] %s" % (return_code, (stderr_text or stdout_text)[-400:]))
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

    performance_values = [run["ns_per_day"] for run in runs]
    mean_value = statistics.mean(performance_values)
    median_value = statistics.median(performance_values)
    std_value = statistics.stdev(performance_values) if len(performance_values) > 1 else 0.0
    cv_value = (std_value / mean_value * 100.0) if mean_value else 0.0
    return dict(
        system=system_name,
        ensemble=ensemble,
        dt=time_step,
        mean=mean_value,
        median=median_value,
        std=std_value,
        cv=cv_value,
        n=len(performance_values),
        runs=runs,
        tuned=tuned_values,
        note=note,
        num_atoms=num_atoms,
        input_options=options,
        autotune_log=autotune_log,
        nbupdate_period=actual_nbupdate_period,
        baroscale_period=actual_baroscale_period,
        mpi_procs=args.mpi_procs,
        omp_threads=args.omp_threads,
    )


def csv_fieldnames(input_columns):
    """Return the CSV header for benchmark result rows."""
    return [
        "system", "ensemble", "dt", "run_id", "ns_per_day", "wall_seconds",
        "ns_per_day_mean", "ns_per_day_median", "ns_per_day_std", "cv_pct",
        "n_measure", "mpi_procs", "omp_threads", "num_atoms",
        "nbupdate_period", "baroscale_period", "tuners", "note",
        "autotune_log", "production_log",
        "tuned_cell_size", "tuned_pairlistdist", "tuned_nbupdate_period",
    ] + ["tuned_%s" % kernel_key for kernel_key in KERNEL_KEYS] + [
        "input_options_json",
    ] + list(input_columns)


def planned_input_columns(system_names, ensembles, time_steps, args):
    """Return input-option CSV columns for the selected benchmark matrix."""
    sentinel_tuned_values = {
        "kernel": dict((kernel_key, 0) for kernel_key in KERNEL_KEYS),
        "cell_size": 0.0,
        "pairlistdist": 0.0,
        "nbupdate_period": 0,
    }
    columns = set()
    for system_name in system_names:
        for ensemble in ensembles:
            for time_step in time_steps:
                tag = "%s_%s_%s" % (system_name, ensemble, time_step)
                base_path = os.path.join(INPUTS, tag + ".inp")
                if not os.path.isfile(base_path):
                    continue
                base_text = open(base_path).read()
                pinned_text = build_pinned_input(
                    base_text, sentinel_tuned_values, args.tune_set, args.num_steps,
                    args.eneout_period, args.nbupdate_period, args.baroscale_period,
                )
                columns.update(input_options(pinned_text))
    return sorted(columns)


def result_csv_rows(result, tune_set, input_columns):
    """Return CSV row dictionaries for one completed benchmark cell."""
    tuned_values = result["tuned"]
    kernel_values = tuned_values.get("kernel") or {}
    options = result["input_options"]
    missing_columns = sorted(set(options) - set(input_columns))
    if missing_columns:
        raise RuntimeError("CSV header is missing input option column(s): %s" %
                           ",".join(missing_columns))

    base_row = {
        "system": result["system"],
        "ensemble": result["ensemble"],
        "dt": result["dt"],
        "ns_per_day_mean": "%.3f" % result["mean"],
        "ns_per_day_median": "%.3f" % result["median"],
        "ns_per_day_std": "%.3f" % result["std"],
        "cv_pct": "%.2f" % result["cv"],
        "n_measure": result["n"],
        "mpi_procs": result["mpi_procs"],
        "omp_threads": result["omp_threads"],
        "num_atoms": result["num_atoms"] or "",
        "nbupdate_period": result["nbupdate_period"] or "",
        "baroscale_period": result["baroscale_period"] or "",
        "tuners": "|".join(sorted(tune_set)) or "none",
        "note": result.get("note", ""),
        "autotune_log": result.get("autotune_log", ""),
        "tuned_cell_size": ("%.3f" % tuned_values["cell_size"]) if "cell_size" in tuned_values else "",
        "tuned_pairlistdist": ("%.3f" % tuned_values["pairlistdist"]) if "pairlistdist" in tuned_values else "",
        "tuned_nbupdate_period": tuned_values.get("nbupdate_period", ""),
        "input_options_json": json.dumps(options, sort_keys=True, separators=(",", ":")),
    }
    for kernel_key in KERNEL_KEYS:
        base_row["tuned_%s" % kernel_key] = kernel_values.get(kernel_key, "")
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


def append_result_csv(csv_path, result, tune_set, input_columns):
    """Append one completed cell to the result CSV."""
    with open(csv_path, "a", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_fieldnames(input_columns), extrasaction="raise")
        writer.writerows(result_csv_rows(result, tune_set, input_columns))
        csv_file.flush()
        os.fsync(csv_file.fileno())


def build_parser():
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--systems", dest="systems_text", metavar="SYSTEMS", default=",".join(ALL_SYSTEMS),
                        help="comma list; default all 8: " + ",".join(ALL_SYSTEMS))
    parser.add_argument("--ensembles", dest="ensembles_text", metavar="ENSEMBLES", default=",".join(ALL_ENSEMBLES),
                        help="nve,nvt,npt")
    parser.add_argument("--dt", dest="time_steps_text", metavar="DT", default="2,4",
                        help="comma list of 2 and/or 4 (fs)")
    parser.add_argument("--warmup", type=int, default=2,
                        help="discarded warm-up runs (default 2)")
    parser.add_argument("--measure", type=int, default=10,
                        help="timed runs (default 10)")
    parser.add_argument("--tune", default="kernel",
                        help="which autotuners: comma list of kernel,cell,nblist | all | none (default kernel)")
    parser.add_argument("--full-autotune", action="store_true",
                        help="shortcut for --tune kernel,cell,nblist")
    parser.add_argument("--nsteps", dest="num_steps", metavar="NSTEPS", type=positive_int, default=100000,
                        help="measurement-run nsteps (default 100000)")
    parser.add_argument("--tune-nsteps", dest="tune_num_steps", metavar="TUNE_NSTEPS",
                        type=positive_int, default=50000,
                        help="tuning-run nsteps (default 50000)")
    parser.add_argument("--eneout-period", type=positive_int, default=1000,
                        help="DYNAMICS eneout_period for tune and measurement inputs (default 1000)")
    parser.add_argument("--mpi-procs", type=positive_int, default=1,
                        help="MPI process count for each spdyn launch (default 1)")
    parser.add_argument("--omp-threads", type=positive_int, default=1,
                        help="OMP_NUM_THREADS for each spdyn launch (default 1)")
    parser.add_argument("--nbupdate-period", type=positive_int, default=None,
                        help="DYNAMICS nbupdate_period override for work inputs")
    parser.add_argument("--baroscale-period", type=positive_int, default=None,
                        help="DYNAMICS baroscale_period override for NPT work inputs")
    parser.add_argument("--timeout", type=int, default=7200,
                        help="per-run timeout seconds (default 7200)")
    parser.add_argument("--lock", default="/tmp/bench.lock",
                        help="advisory serialisation lockfile")
    parser.add_argument("--allow-failures", action="store_true",
                        help="write a partial CSV and exit 0 even if selected cells fail")
    parser.add_argument("--timestamp", default=None,
                        help="results filename stamp (else auto)")
    parser.add_argument("--out", default=None,
                        help="explicit CSV output path (overrides --timestamp)")
    parser.add_argument("--spdyn", default=SPDYN,
                        help="path to the spdyn binary")
    return parser


def resolve_tune_set(args, parser):
    """Return the selected autotuner set."""
    if args.full_autotune:
        return {"kernel", "cell", "nblist"}

    raw_tuners = [tuner.strip().lower() for tuner in args.tune.split(",") if tuner.strip()]
    if "all" in raw_tuners:
        return {"kernel", "cell", "nblist"}
    if "none" in raw_tuners or not raw_tuners:
        return set()

    aliases = {"nbupdate": "nblist", "nb": "nblist", "cell_size": "cell"}
    tune_set = set(aliases.get(tuner, tuner) for tuner in raw_tuners)
    unknown_tuners = tune_set - {"kernel", "cell", "nblist"}
    if unknown_tuners:
        parser.error("unknown tuner(s): %s (use kernel,cell,nblist,all,none)" %
                     ",".join(sorted(unknown_tuners)))
    return tune_set


def parse_time_steps(time_steps_text, parser):
    """Return selected time-step labels from the --dt option."""
    time_steps = []
    for time_step_text in time_steps_text.split(","):
        normalized_time_step = time_step_text.strip().replace("fs", "")
        if normalized_time_step in ("2", "4"):
            time_steps.append(normalized_time_step + "fs")
    if not time_steps:
        parser.error("--dt must include 2 and/or 4")
    return time_steps


def parse_systems(systems_text, parser):
    """Return selected system names from the --systems option."""
    system_names = [system_name.strip() for system_name in systems_text.split(",") if system_name.strip()]
    for system_name in system_names:
        if system_name not in ALL_SYSTEMS:
            parser.error("unknown system %r (choices: %s)" % (system_name, ",".join(ALL_SYSTEMS)))
    return system_names


def parse_ensembles(ensembles_text):
    """Return selected ensemble names from the --ensembles option."""
    return [ensemble.strip().lower() for ensemble in ensembles_text.split(",") if ensemble.strip()]


def main():
    """Run the benchmark command-line program."""
    global SPDYN

    parser = build_parser()
    args = parser.parse_args()
    validate_eneout_period(args, parser)

    SPDYN = os.path.abspath(os.path.expanduser(args.spdyn))
    if not (os.path.isfile(SPDYN) and os.access(SPDYN, os.X_OK)):
        sys.exit("error: spdyn not found or not executable: %s" % SPDYN)

    tune_set = resolve_tune_set(args, parser)
    args.tune_set = tune_set
    if "nblist" in tune_set and args.nbupdate_period is not None:
        parser.error("--nbupdate-period cannot be combined with nblist autotune; "
                     "pairlistdist and nbupdate_period are tuned as a coupled pair")

    system_names = parse_systems(args.systems_text, parser)
    ensembles = parse_ensembles(args.ensembles_text)
    time_steps = parse_time_steps(args.time_steps_text, parser)

    results = []
    failures = []
    csv_path = None
    log_root = None
    log("waiting for benchmark lock: %s" % args.lock)
    lock_fd = acquire_lock(args.lock, "run_benchmark")
    try:
        try:
            extracted_systems = ensure_system_data(system_names)
            timestamp = args.timestamp or datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
            requested_csv_path = args.out or os.path.join(RESULTS, timestamp + ".csv")
            csv_path, log_root = unique_run_paths(requested_csv_path)
            args.log_root = log_root
            args.input_log_dir = os.path.join(log_root, "inputs")
            os.makedirs(os.path.join(log_root, "autotune"), exist_ok=True)
            os.makedirs(os.path.join(log_root, "production"), exist_ok=True)
            os.makedirs(args.input_log_dir, exist_ok=True)
            open_benchmark_log(os.path.join(log_root, "benchmark.log"))

            input_columns = planned_input_columns(system_names, ensembles, time_steps, args)
            create_results_csv(csv_path, input_columns)

            log("=== GENESIS GPU benchmark ===")
            log("spdyn      : %s" % SPDYN)
            log("systems    : %s" % ",".join(system_names))
            log("ensembles  : %s" % ",".join(ensembles))
            log("dt         : %s" % ",".join(time_steps))
            log("tuners     : %s" % (",".join(sorted(tune_set)) or "none"))
            log("warmup/measure: %d/%d   nsteps meas/tune: %d/%d" %
                (args.warmup, args.measure, args.num_steps, args.tune_num_steps))
            log("eneout     : %d" % args.eneout_period)
            log("mpi/omp    : %d/%d" % (args.mpi_procs, args.omp_threads))
            log("nbupdate  : %s" % (args.nbupdate_period if args.nbupdate_period is not None else "input"))
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

            for system_name in system_names:
                for ensemble in ensembles:
                    for time_step in time_steps:
                        result = run_cell(system_name, ensemble, time_step, args)
                        if result:
                            results.append(result)
                            append_result_csv(csv_path, result, tune_set, input_columns)
                        else:
                            failures.append("%s_%s_%s" % (system_name, ensemble, time_step))
        finally:
            release_lock(lock_fd)

        table_lines = summary_table_lines(results, tune_set, args.measure, args.num_steps)
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
