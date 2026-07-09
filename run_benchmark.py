#!/usr/bin/env python3
"""GPU performance benchmark driver for GENESIS spdyn.

Protocol, per selected (system, ensemble, dt) cell:
  1. TUNE   : run ONCE with the requested autotuners ON, parse the [AUTOTUNE] report
              to extract the tuned kernel block sizes / cell_size / pairlistdist / nbupdate_period.
              A manual --nbupdate-period is incompatible with the nblist tuner because
              pairlistdist and nbupdate_period are selected as a coupled pair.
  2. PIN    : write a FRESH input that hard-codes the tuned values and turns ALL autotuners OFF,
              so the measured runs carry no autotune instrumentation.
  3. WARM-UP: run the pinned input N times (default 2), discard timings (heat the GPU / stabilise clocks).
  4. MEASURE: run the pinned input M times (default 10), collect the [PERFORMANCE] ns/day each.
  5. Report mean, median, std, and cv ns/day.

Every spdyn launch is `mpirun -np <mpi-procs> spdyn <inp>` from the benchmark root with
  GENESIS_GPU_PROFILE=0 OMP_NUM_THREADS=<omp-threads> HWLOC_COMPONENTS=x86
and the profiler-independent ns/day is read from the [PERFORMANCE] line.

Input data are stored as data/<system>.tgz archives in git. If data/<system>/ is
missing, the selected system archive is extracted while the benchmark lock is held.

Runs are serialised via an advisory lockfile (default /tmp/bench.lock) so concurrent
agents never benchmark at the same time on the shared machine. The file may remain
on disk; only a live flock owner blocks another run.

Examples:
  # default tune (kernel only), 2 systems, small window, quick:
  python3 run_benchmark.py --systems dhfr,apoa1 --ensembles nve,nvt --dt 2 --warmup 1 --measure 3 --nsteps 2000 --tune-nsteps 3000
  # override neighbor-list/barostat scaling periods for generated work inputs:
  python3 run_benchmark.py --systems dhfr --nbupdate-period 20 --baroscale-period 10
  # choose MPI/OpenMP parallelism for every spdyn launch:
  python3 run_benchmark.py --systems dhfr --mpi-procs 1 --omp-threads 1
  # full run, all 8 systems x 3 ensembles x 2 dt, all three autotuners:
  python3 run_benchmark.py --full-autotune
  # choose a subset of tuners:
  python3 run_benchmark.py --tune kernel,nblist
See README.md for the full protocol.
"""
import os, re, sys, time, argparse, statistics, subprocess, datetime, fcntl, signal, tarfile, hashlib, shutil, tempfile, csv, shlex, json

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
INPUTS = os.path.join(HERE, "inputs")
RESULTS = os.path.join(HERE, "results")
SPDYN = os.path.normpath(os.path.join(HERE, "..", "src", "spdyn_singlempi", "spdyn"))

ALL_SYSTEMS = ["dhfr", "apoa1", "uun", "factorix", "bpti", "dppc", "ake", "stmv"]
ALL_ENSEMBLES = ["nve", "nvt", "npt"]
ALL_DT = ["2fs", "4fs"]
ARCHIVE_SENTINEL = ".archive_sha256"
KERNEL_KEYS = ["kernel_pme_spread_block", "kernel_pme_influence_block", "kernel_bonded_block",
               "kernel_constraints_block", "kernel_nonbond_inter_minblocks", "kernel_nonbond_intra_minblocks"]

PERF_RE = re.compile(r'\[PERFORMANCE\]\s*(?:performance:\s*)?([0-9.]+)\s*ns/day')
NUM_ATOMS_RES = [
    re.compile(r'^\s*num_atoms\s*=\s*([0-9]+)\b', re.IGNORECASE | re.MULTILINE),
    re.compile(r'^\s*natom\s*=\s*([0-9]+)\b', re.IGNORECASE | re.MULTILINE),
    re.compile(r'^\s*number\s+of\s+atoms\s*[:=]\s*([0-9]+)\b', re.IGNORECASE | re.MULTILINE),
]

# ---------------------------------------------------------------------------
# GENESIS control-file editing (section-aware; keys/sections are case-insensitive)
# ---------------------------------------------------------------------------
SECTION_RE = re.compile(r'^\s*\[([A-Za-z0-9_]+)\]\s*$')
KEYLINE_RE = re.compile(r'^\s*([A-Za-z_]\w*)\s*=')


def split_sections(text):
    """-> list of [section_name_or_None, [raw lines incl. header]]."""
    blocks, name, cur = [], None, []
    for ln in text.splitlines():
        m = SECTION_RE.match(ln)
        if m:
            blocks.append([name, cur])
            name, cur = m.group(1), [ln]
        else:
            cur.append(ln)
    blocks.append([name, cur])
    return blocks


def join_sections(blocks):
    out = []
    for _, lines in blocks:
        out.extend(lines)
    return "\n".join(out).rstrip("\n") + "\n"


def _find_block(blocks, section):
    for b in blocks:
        if b[0] and b[0].upper() == section.upper():
            return b
    return None


def set_kv(blocks, section, key, value):
    """Set key=value inside [section]; create the section/line if absent. Continuation- and
    comment-safe (skips '\\'-continued and '#'-commented lines)."""
    b = _find_block(blocks, section)
    val = str(value)
    if b is None:
        blocks.append([section, ["[%s]" % section, "%-16s = %s" % (key, val)]])
        return
    lines = b[1]
    prev_cont = False
    for i in range(1, len(lines)):
        ln = lines[i]
        stripped = ln.strip()
        is_cont = prev_cont
        prev_cont = stripped.endswith("\\")
        if is_cont or stripped.startswith("#"):
            continue
        m = KEYLINE_RE.match(ln)
        if m and m.group(1).lower() == key.lower():
            lines[i] = "%-16s = %s" % (key, val)
            return
    # not found: insert right after the header
    lines.insert(1, "%-16s = %s" % (key, val))


def has_kv(blocks, section, key):
    b = _find_block(blocks, section)
    if b is None:
        return False
    prev_cont = False
    for ln in b[1][1:]:
        stripped = ln.strip()
        is_cont = prev_cont
        prev_cont = stripped.endswith("\\")
        if is_cont or stripped.startswith("#"):
            continue
        m = KEYLINE_RE.match(ln)
        if m and m.group(1).lower() == key.lower():
            return True
    return False


def get_kv(text, section, key):
    b = _find_block(split_sections(text), section)
    if b is None:
        return None
    prev_cont = False
    for ln in b[1][1:]:
        stripped = ln.strip()
        is_cont = prev_cont
        prev_cont = stripped.endswith("\\")
        if is_cont or stripped.startswith("#"):
            continue
        m = KEYLINE_RE.match(ln)
        if m and m.group(1).lower() == key.lower():
            return ln.split("=", 1)[1].split("#", 1)[0].strip()
    return None


def set_window(blocks, nsteps, eneout_period):
    """Set nsteps and energy-output period for a run window."""
    set_kv(blocks, "DYNAMICS", "nsteps", nsteps)
    set_kv(blocks, "DYNAMICS", "eneout_period", eneout_period)


def set_period_overrides(blocks, nbupdate_period, baroscale_period):
    if nbupdate_period is not None:
        set_kv(blocks, "DYNAMICS", "nbupdate_period", nbupdate_period)
    if baroscale_period is not None and has_kv(blocks, "DYNAMICS", "barostat_period"):
        set_kv(blocks, "DYNAMICS", "baroscale_period", baroscale_period)


# ---------------------------------------------------------------------------
# [AUTOTUNE] report parsing
# ---------------------------------------------------------------------------
def parse_autotune(text):
    """Extract tuned values from the [AUTOTUNE] report. Returns a dict with any of:
    'kernel' (dict), 'cell_size' (float), 'pairlistdist' (float), 'nbupdate_period' (int).
    Missing sections are simply absent."""
    res = {}
    lines = text.splitlines()

    # Neighborlist: candidate blocks; find the one whose Time line says "(selected)".
    # Read only THIS block's pairlistdist/nbupdate_period -- stop at the next candidate
    # or the start of the cell/kernel section, so we don't read the next block's values.
    for i, ln in enumerate(lines):
        if "ms/window" in ln and "(selected)" in ln:
            pld = nbu = None
            for j in range(i + 1, len(lines)):
                lj = lines[j]
                if "ms/window" in lj or "Cell size autotune" in lj or "Selected" in lj:
                    break
                if pld is None:
                    m = re.search(r'pairlistdist\s*=\s*([0-9.]+)', lj)
                    if m:
                        pld = float(m.group(1))
                if nbu is None:
                    m = re.search(r'nbupdate_period\s*=\s*([0-9]+)', lj)
                    if m:
                        nbu = int(m.group(1))
                if pld is not None and nbu is not None:
                    break
            if pld is not None:
                res["pairlistdist"] = pld
            if nbu is not None:
                res["nbupdate_period"] = nbu
            break

    # Cell size: the line after "Selected configuration:".
    for i, ln in enumerate(lines):
        if "Selected configuration:" in ln:
            for j in range(i + 1, min(i + 4, len(lines))):
                m = re.search(r'cell_size\s*=\s*([0-9.]+)', lines[j])
                if m:
                    res["cell_size"] = float(m.group(1))
                    break
            break

    # Kernel block sizes.
    kern = {}
    for ln in lines:
        m = re.match(r'\s*(kernel_\w+)\s*=\s*([0-9]+)', ln)
        if m and m.group(1) in KERNEL_KEYS:
            kern[m.group(1)] = int(m.group(2))
    if kern:
        res["kernel"] = kern
    return res


def parse_perf(text):
    m = PERF_RE.search(text)
    return float(m.group(1)) if m else None


def parse_num_atoms(text):
    for regex in NUM_ATOMS_RES:
        m = regex.search(text)
        if m:
            return int(m.group(1))
    return None


def _normalize_csv_key(text):
    return re.sub(r'[^0-9a-zA-Z_]+', '_', text.strip().lower()).strip("_")


def input_options(text):
    """Return section-qualified input parameters from a GENESIS control file."""
    opts = {}
    for section, lines in split_sections(text):
        if not section:
            continue
        sec = _normalize_csv_key(section)
        i = 1
        while i < len(lines):
            ln = lines[i]
            stripped = ln.strip()
            if not stripped or stripped.startswith("#"):
                i += 1
                continue
            m = KEYLINE_RE.match(ln)
            if not m:
                i += 1
                continue
            key = _normalize_csv_key(m.group(1))
            value_parts = []
            raw_value = ln.split("=", 1)[1]
            while True:
                raw_no_comment = raw_value.split("#", 1)[0].strip()
                continued = raw_no_comment.endswith("\\")
                if continued:
                    raw_no_comment = raw_no_comment[:-1].rstrip()
                if raw_no_comment:
                    value_parts.append(raw_no_comment)
                i += 1
                if not continued or i >= len(lines):
                    break
                raw_value = lines[i].strip()
            opts["input_%s_%s" % (sec, key)] = " ".join(value_parts)
    return opts


def validate_tuned(tuned, tune_set):
    """Return a list of requested autotune values missing from the parsed report."""
    missing = []
    if "kernel" in tune_set:
        kern = tuned.get("kernel") or {}
        missing.extend(k for k in KERNEL_KEYS if k not in kern)
    if "cell" in tune_set and "cell_size" not in tuned:
        missing.append("cell_size")
    if "nblist" in tune_set:
        if "pairlistdist" not in tuned:
            missing.append("pairlistdist")
        if "nbupdate_period" not in tuned:
            missing.append("nbupdate_period")
    return missing


# ---------------------------------------------------------------------------
# Running spdyn (lock-serialised)
# ---------------------------------------------------------------------------
def acquire_lock(lockpath, label, poll=2.0):
    """Acquire an advisory benchmark lock.

    The old driver used O_EXCL create/remove. If the Python process was killed,
    the file stayed behind and blocked all future runs. flock ties ownership to
    this open fd instead, so a dead process cannot leave a live lock.
    """
    lockdir = os.path.dirname(os.path.abspath(lockpath))
    if lockdir:
        os.makedirs(lockdir, exist_ok=True)
    fd = os.open(lockpath, os.O_CREAT | os.O_RDWR, 0o666)
    while True:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            os.ftruncate(fd, 0)
            os.lseek(fd, 0, os.SEEK_SET)
            os.write(fd, ("%d %s\n" % (os.getpid(), label)).encode())
            return fd
        except BlockingIOError:
            time.sleep(poll)


def release_lock(lock_fd):
    try:
        os.ftruncate(lock_fd, 0)
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
    finally:
        os.close(lock_fd)


def _safe_tar_member(data_root, sysname, member):
    """Reject archive members that could write outside data/<system>/."""
    name = member.name
    expected = sysname + "/"
    if name == sysname:
        if not member.isdir():
            return False
    elif not name.startswith(expected):
        return False
    target = os.path.abspath(os.path.join(data_root, name))
    system_root = os.path.abspath(os.path.join(data_root, sysname))
    if target != system_root and not target.startswith(system_root + os.sep):
        return False
    # The benchmark archives are plain files/directories. Reject links and
    # special files so extraction cannot depend on external filesystem state.
    if not (member.isfile() or member.isdir()):
        return False
    return True


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def system_data_is_current(system_dir, regular_members, archive_digest):
    sentinel = os.path.join(system_dir, ARCHIVE_SENTINEL)
    try:
        if open(sentinel).read().strip() != archive_digest:
            return False
    except OSError:
        return False
    for member in regular_members:
        path = os.path.join(DATA, member.name)
        if not os.path.isfile(path):
            return False
        if os.path.getsize(path) != member.size:
            return False
    return True


def extract_system_archive(archive, sysname, members, archive_digest):
    tmp_root = tempfile.mkdtemp(prefix=".extract-%s-" % sysname, dir=DATA)
    tmp_system_dir = os.path.join(tmp_root, sysname)
    system_dir = os.path.join(DATA, sysname)
    try:
        with tarfile.open(archive, "r:gz") as tf:
            tf.extractall(tmp_root, members)
        if not os.path.isdir(tmp_system_dir):
            raise RuntimeError("%s did not contain %s/" % (archive, sysname))
        with open(os.path.join(tmp_system_dir, ARCHIVE_SENTINEL), "w") as f:
            f.write(archive_digest + "\n")
        if os.path.exists(system_dir):
            shutil.rmtree(system_dir)
        os.rename(tmp_system_dir, system_dir)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)


def ensure_system_data(systems):
    """Extract data/<system>.tgz for selected systems whose data dir is absent."""
    extracted = []
    os.makedirs(DATA, exist_ok=True)
    for sysname in systems:
        system_dir = os.path.join(DATA, sysname)
        archive = os.path.join(DATA, sysname + ".tgz")
        if not os.path.isfile(archive):
            raise RuntimeError("missing %s and %s" % (system_dir, archive))
        archive_digest = sha256_file(archive)
        with tarfile.open(archive, "r:gz") as tf:
            members = tf.getmembers()
            bad = [m.name for m in members if not _safe_tar_member(DATA, sysname, m)]
            if bad:
                raise RuntimeError("unsafe member(s) in %s: %s" %
                                   (archive, ", ".join(bad[:5])))
        regular_members = [m for m in members if m.isfile()]
        if os.path.isdir(system_dir) and system_data_is_current(system_dir, regular_members, archive_digest):
            continue
        extract_system_archive(archive, sysname, members, archive_digest)
        extracted.append(sysname)
    return extracted


def terminate_process_group(proc):
    if proc is None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5.0)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()


def write_spdyn_log(log_path, phase, tag, run_id, inp_path, cmd, env, rc, wall, timeout, out, err):
    if log_path is None:
        return
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as f:
        f.write("# GENESIS benchmark run log\n")
        f.write("# phase: %s\n" % (phase or "unknown"))
        f.write("# tag: %s\n" % (tag or "unknown"))
        f.write("# run_id: %s\n" % (run_id if run_id is not None else ""))
        f.write("# input: %s\n" % inp_path)
        f.write("# cwd: %s\n" % HERE)
        f.write("# command: %s\n" % " ".join(shlex.quote(str(x)) for x in cmd))
        f.write("# env: GENESIS_GPU_PROFILE=%s OMP_NUM_THREADS=%s HWLOC_COMPONENTS=%s\n" %
                (env.get("GENESIS_GPU_PROFILE", ""), env.get("OMP_NUM_THREADS", ""),
                 env.get("HWLOC_COMPONENTS", "")))
        f.write("# return_code: %s\n" % rc)
        f.write("# wall_seconds: %.6f\n" % wall)
        f.write("# timeout_seconds: %s\n" % timeout)
        f.write("\n# stdout\n")
        f.write(out or "")
        if out and not out.endswith("\n"):
            f.write("\n")
        f.write("\n# stderr\n")
        f.write(err or "")
        if err and not err.endswith("\n"):
            f.write("\n")


def run_spdyn(inp_path, timeout, mpi_procs, omp_threads, log_path=None,
              phase=None, tag=None, run_id=None):
    env = dict(os.environ, GENESIS_GPU_PROFILE="0", OMP_NUM_THREADS=str(omp_threads),
               HWLOC_COMPONENTS="x86")
    cmd = ["mpirun", "-np", str(mpi_procs), SPDYN, inp_path]
    p = None
    out = ""
    err = ""
    rc = None
    t0 = time.time()
    try:
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, env=env, cwd=HERE, start_new_session=True)
        try:
            out, err = p.communicate(timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            terminate_process_group(p)
            out, err = p.communicate()
            rc = 124
        wall = time.time() - t0
    except KeyboardInterrupt:
        terminate_process_group(p)
        if p is not None:
            try:
                out, err = p.communicate(timeout=1.0)
            except subprocess.TimeoutExpired:
                terminate_process_group(p)
                out, err = p.communicate()
            except BaseException:
                out = out or ""
                err = err or ""
        rc = 130
        wall = time.time() - t0
        err = (err or "") + "\n[run_benchmark] interrupted by KeyboardInterrupt\n"
        write_spdyn_log(log_path, phase, tag, run_id, inp_path, cmd, env, rc, wall, timeout, out, err)
        raise
    write_spdyn_log(log_path, phase, tag, run_id, inp_path, cmd, env, rc, wall, timeout, out, err)
    return out, err, rc, wall


# ---------------------------------------------------------------------------
# Tune / pin input builders
# ---------------------------------------------------------------------------
def build_tune_input(base_text, tune_set, nsteps, eneout_period, nbupdate_period, baroscale_period):
    blocks = split_sections(base_text)
    set_window(blocks, nsteps, eneout_period)
    set_period_overrides(blocks, nbupdate_period, baroscale_period)
    if "cell" in tune_set:
        set_kv(blocks, "ENERGY", "cell_size_autotune", "YES")
    if "nblist" in tune_set:
        set_kv(blocks, "DYNAMICS", "nbupdate_autotune", "YES")
    if "kernel" in tune_set:
        set_kv(blocks, "GPU", "kernel_autotune", "YES")
    return join_sections(blocks)


def build_pinned_input(base_text, tuned, tune_set, nsteps, eneout_period, nbupdate_period, baroscale_period):
    blocks = split_sections(base_text)
    set_window(blocks, nsteps, eneout_period)
    # All autotuners explicitly OFF.
    set_kv(blocks, "ENERGY", "cell_size_autotune", "NO")
    set_kv(blocks, "DYNAMICS", "nbupdate_autotune", "NO")
    set_kv(blocks, "GPU", "kernel_autotune", "NO")
    # Pin whatever was requested AND parsed.
    if "kernel" in tune_set and tuned.get("kernel"):
        for k in KERNEL_KEYS:
            if k in tuned["kernel"]:
                set_kv(blocks, "GPU", k, tuned["kernel"][k])
    if "cell" in tune_set and "cell_size" in tuned:
        set_kv(blocks, "ENERGY", "cell_size", "%.3f" % tuned["cell_size"])
    if "nblist" in tune_set:
        if "pairlistdist" in tuned:
            set_kv(blocks, "ENERGY", "pairlistdist", "%.3f" % tuned["pairlistdist"])
        if "nbupdate_period" in tuned:
            set_kv(blocks, "DYNAMICS", "nbupdate_period", tuned["nbupdate_period"])
    set_period_overrides(blocks, nbupdate_period, baroscale_period)
    return join_sections(blocks)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
BENCHMARK_FH = None
BENCHMARK_BUFFER = []


def open_benchmark_log(path):
    global BENCHMARK_FH, BENCHMARK_BUFFER
    os.makedirs(os.path.dirname(path), exist_ok=True)
    BENCHMARK_FH = open(path, "w")
    for msg in BENCHMARK_BUFFER:
        BENCHMARK_FH.write(msg + "\n")
    BENCHMARK_FH.flush()
    BENCHMARK_BUFFER = []


def close_benchmark_log():
    global BENCHMARK_FH
    if BENCHMARK_FH is not None:
        BENCHMARK_FH.close()
        BENCHMARK_FH = None


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()
    if BENCHMARK_FH is not None:
        BENCHMARK_FH.write(msg + "\n")
        BENCHMARK_FH.flush()
    else:
        BENCHMARK_BUFFER.append(msg)


def emit(msg=""):
    print(msg)
    if BENCHMARK_FH is not None:
        BENCHMARK_FH.write(msg + "\n")
        BENCHMARK_FH.flush()


def write_summary_log(path, lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")


def summary_table_lines(results, tune, measure, nsteps):
    lines = [
        "=== ns/day (mean/median +- std, cv%%) : tuners=%s, measure=%d, nsteps=%d ===" %
        (",".join(sorted(tune)) or "none", measure, nsteps),
        "%-10s %-5s %-4s %12s %12s %10s %7s  %s" %
        ("system", "ens", "dt", "mean", "median", "+-std", "cv%", "note"),
        "-" * 82,
    ]
    for r in results:
        lines.append("%-10s %-5s %-4s %12.2f %12.2f %10.2f %6.1f%%  %s" %
                     (r["system"], r["ensemble"], r["dt"], r["mean"], r["median"],
                      r["std"], r["cv"], r.get("note", "")))
    return lines


def unique_run_paths(path):
    """Pick non-existing CSV and log-directory paths while the suite lock is held."""
    root, ext = os.path.splitext(path)
    log_root = root if ext.lower() == ".csv" else path + ".logs"
    if not os.path.exists(path) and not os.path.exists(log_root):
        return path, log_root
    for i in range(1, 1000):
        cand_csv = "%s-%d%s" % (root, i, ext)
        cand_log_root = "%s-%d" % (log_root, i)
        if not os.path.exists(cand_csv) and not os.path.exists(cand_log_root):
            return cand_csv, cand_log_root
    raise RuntimeError("could not find unused output/log paths for %s" % path)


def run_log_path(log_root, subdir, tag, run_name):
    return os.path.join(log_root, subdir, "%s_%s.log" % (tag, run_name))


def positive_int(text):
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if value < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return value


def validate_eneout_period(args, parser):
    if args.nsteps % args.eneout_period != 0:
        parser.error("--nsteps (%d) must be a multiple of --eneout-period (%d)" %
                     (args.nsteps, args.eneout_period))
    if args.tune_nsteps % args.eneout_period != 0:
        parser.error("--tune-nsteps (%d) must be a multiple of --eneout-period (%d)" %
                     (args.tune_nsteps, args.eneout_period))


def run_cell(sysname, ens, dt, args):
    tag = "%s_%s_%s" % (sysname, ens, dt)
    base_path = os.path.join(INPUTS, tag + ".inp")
    if not os.path.isfile(base_path):
        log("  [skip] no base input %s" % base_path)
        return None
    base_text = open(base_path).read()
    tuned = {}
    note = ""
    num_atoms = None
    autotune_log = ""

    # 1. TUNE
    if args.tune_set:
        tune_inp = os.path.join(args.input_log_dir, tag + ".tune.inp")
        open(tune_inp, "w").write(build_tune_input(base_text, args.tune_set, args.tune_nsteps,
                                                   args.eneout_period,
                                                   args.nbupdate_period, args.baroscale_period))
        tune_log = run_log_path(args.log_root, "autotune", tag, 1)
        autotune_log = os.path.relpath(tune_log, HERE)
        log("  [tune] %s  tuners=%s  nsteps=%d" % (tag, ",".join(sorted(args.tune_set)), args.tune_nsteps))
        out, err, rc, wall = run_spdyn(tune_inp, args.timeout, args.mpi_procs, args.omp_threads,
                                       tune_log, "autotune", tag, 1)
        num_atoms = parse_num_atoms((out or "") + "\n" + (err or "")) or num_atoms
        if rc != 0:
            # Only the known cell/nblist runtime re-decomposition overflow is recoverable.
            # Kernel-only tuning failures mean the requested protocol did not run.
            reason = "cell overflow" if "cell too large" in (out + err) else ("rc=%d" % rc)
            if args.tune_set <= {"kernel"} or "cell overflow" not in reason:
                log("  [tune FAILED: %s] %s" % (reason, (err or out)[-400:]))
                return None
            log("  [tune FAILED: %s -> measuring with DEFAULTS] %s" % (reason, tag))
            note = "tune-failed(%s)" % reason
            tuned = {}
        else:
            tuned = parse_autotune(out)
            missing = validate_tuned(tuned, args.tune_set)
            if missing:
                log("  [tune FAILED: missing parsed value(s): %s] %s" %
                    (",".join(missing), tag))
                return None
            log("    tuned: " + (", ".join("%s=%s" % (k, v) for k, v in tuned.items() if k != "kernel")
                                 + ((" kernel=" + str(tuned["kernel"])) if "kernel" in tuned else "")
                                 or "(nothing parsed)"))

    # 2. PIN
    pin_inp = os.path.join(args.input_log_dir, tag + ".pinned.inp")
    pin_text = build_pinned_input(base_text, tuned, args.tune_set, args.nsteps,
                                  args.eneout_period,
                                  args.nbupdate_period, args.baroscale_period)
    open(pin_inp, "w").write(pin_text)
    actual_nbupdate_period = get_kv(pin_text, "DYNAMICS", "nbupdate_period")
    actual_baroscale_period = get_kv(pin_text, "DYNAMICS", "baroscale_period")
    opts = input_options(pin_text)

    # 3. WARM-UP
    for w in range(args.warmup):
        log("  [warmup %d/%d] %s" % (w + 1, args.warmup, tag))
        warmup_log = run_log_path(args.log_root, "production", tag, "warmup_%d" % (w + 1))
        out, err, rc, wall = run_spdyn(pin_inp, args.timeout, args.mpi_procs, args.omp_threads,
                                       warmup_log, "production-warmup", tag, w + 1)
        num_atoms = parse_num_atoms((out or "") + "\n" + (err or "")) or num_atoms
        if rc != 0:
            log("  [warmup FAILED rc=%d] %s" % (rc, (err or out)[-400:]))
            return None

    # 4. MEASURE
    runs = []
    for m in range(args.measure):
        production_log = run_log_path(args.log_root, "production", tag, m + 1)
        out, err, rc, wall = run_spdyn(pin_inp, args.timeout, args.mpi_procs, args.omp_threads,
                                       production_log, "production-measure", tag, m + 1)
        num_atoms = parse_num_atoms((out or "") + "\n" + (err or "")) or num_atoms
        if rc != 0:
            log("  [measure FAILED rc=%d] %s" % (rc, (err or out)[-400:]))
            return None
        v = parse_perf(out)
        if v is None:
            log("  [measure FAILED: no [PERFORMANCE] parsed] %s" % ((err or out)[-400:]))
            return None
        runs.append(dict(run_id=m + 1, ns_per_day=v, wall_seconds=wall,
                         production_log=os.path.relpath(production_log, HERE)))
        log("  [measure %d/%d] %s  ns/day=%s  (wall %.1fs)" % (m + 1, args.measure, tag,
                                                               ("%.2f" % v) if v else "None", wall))
    vals = [r["ns_per_day"] for r in runs]
    mean = statistics.mean(vals)
    median = statistics.median(vals)
    std = statistics.stdev(vals) if len(vals) > 1 else 0.0
    cv = (std / mean * 100.0) if mean else 0.0
    return dict(system=sysname, ensemble=ens, dt=dt, mean=mean, median=median, std=std, cv=cv,
                n=len(vals), runs=runs, tuned=tuned, note=note, num_atoms=num_atoms,
                input_options=opts, autotune_log=autotune_log,
                nbupdate_period=actual_nbupdate_period,
                baroscale_period=actual_baroscale_period,
                mpi_procs=args.mpi_procs, omp_threads=args.omp_threads)


def csv_fieldnames(input_columns):
    return [
        "system", "ensemble", "dt", "run_id", "ns_per_day", "wall_seconds",
        "ns_per_day_mean", "ns_per_day_median", "ns_per_day_std", "cv_pct",
        "n_measure", "mpi_procs", "omp_threads", "num_atoms",
        "nbupdate_period", "baroscale_period", "tuners", "note",
        "autotune_log", "production_log",
        "tuned_cell_size", "tuned_pairlistdist", "tuned_nbupdate_period",
    ] + ["tuned_%s" % k for k in KERNEL_KEYS] + [
        "input_options_json",
    ] + list(input_columns)


def planned_input_columns(systems, ensembles, dts, args):
    sentinel = {
        "kernel": dict((k, 0) for k in KERNEL_KEYS),
        "cell_size": 0.0,
        "pairlistdist": 0.0,
        "nbupdate_period": 0,
    }
    columns = set()
    for sysname in systems:
        for ens in ensembles:
            for dt in dts:
                tag = "%s_%s_%s" % (sysname, ens, dt)
                base_path = os.path.join(INPUTS, tag + ".inp")
                if not os.path.isfile(base_path):
                    continue
                base_text = open(base_path).read()
                pin_text = build_pinned_input(base_text, sentinel, args.tune_set, args.nsteps,
                                              args.eneout_period,
                                              args.nbupdate_period, args.baroscale_period)
                columns.update(input_options(pin_text))
    return sorted(columns)


def result_csv_rows(result, tune, input_columns):
    t = result["tuned"]
    kern = t.get("kernel") or {}
    options = result["input_options"]
    missing = sorted(set(options) - set(input_columns))
    if missing:
        raise RuntimeError("CSV header is missing input option column(s): %s" %
                           ",".join(missing))
    base = {
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
        "tuners": "|".join(sorted(tune)) or "none",
        "note": result.get("note", ""),
        "autotune_log": result.get("autotune_log", ""),
        "tuned_cell_size": ("%.3f" % t["cell_size"]) if "cell_size" in t else "",
        "tuned_pairlistdist": ("%.3f" % t["pairlistdist"]) if "pairlistdist" in t else "",
        "tuned_nbupdate_period": t.get("nbupdate_period", ""),
        "input_options_json": json.dumps(options, sort_keys=True, separators=(",", ":")),
    }
    for k in KERNEL_KEYS:
        base["tuned_%s" % k] = kern.get(k, "")
    base.update(options)
    rows = []
    for run in result["runs"]:
        row = dict(base)
        row.update({
            "run_id": run["run_id"],
            "ns_per_day": "%.3f" % run["ns_per_day"],
            "wall_seconds": "%.3f" % run["wall_seconds"],
            "production_log": run["production_log"],
        })
        rows.append(row)
    return rows


def create_results_csv(csv_path, input_columns):
    with open(csv_path, "x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fieldnames(input_columns))
        writer.writeheader()


def append_result_csv(csv_path, result, tune, input_columns):
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fieldnames(input_columns), extrasaction="raise")
        writer.writerows(result_csv_rows(result, tune, input_columns))
        f.flush()
        os.fsync(f.fileno())


def write_results_csv(csv_path, results, tune):
    input_columns = sorted({k for r in results for k in r["input_options"]})
    columns = csv_fieldnames(input_columns)
    with open(csv_path, "x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerows(result_csv_rows(r, tune, input_columns))


def main():
    global SPDYN
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--systems", default=",".join(ALL_SYSTEMS),
                    help="comma list; default all 8: " + ",".join(ALL_SYSTEMS))
    ap.add_argument("--ensembles", default=",".join(ALL_ENSEMBLES), help="nve,nvt,npt")
    ap.add_argument("--dt", default="2,4", help="comma list of 2 and/or 4 (fs)")
    ap.add_argument("--warmup", type=int, default=2, help="discarded warm-up runs (default 2)")
    ap.add_argument("--measure", type=int, default=10, help="timed runs (default 10)")
    ap.add_argument("--tune", default="kernel",
                    help="which autotuners: comma list of kernel,cell,nblist | all | none (default kernel)")
    ap.add_argument("--full-autotune", action="store_true",
                    help="shortcut for --tune kernel,cell,nblist")
    ap.add_argument("--nsteps", type=positive_int, default=100000,
                    help="measurement-run nsteps (default 100000)")
    ap.add_argument("--tune-nsteps", type=positive_int, default=50000,
                    help="tuning-run nsteps (default 50000)")
    ap.add_argument("--eneout-period", type=positive_int, default=1000,
                    help="DYNAMICS eneout_period for tune and measurement inputs (default 1000)")
    ap.add_argument("--mpi-procs", type=positive_int, default=1,
                    help="MPI process count for each spdyn launch (default 1)")
    ap.add_argument("--omp-threads", type=positive_int, default=1,
                    help="OMP_NUM_THREADS for each spdyn launch (default 1)")
    ap.add_argument("--nbupdate-period", type=positive_int, default=None,
                    help="DYNAMICS nbupdate_period override for work inputs (default: input value, 10 in generated inputs)")
    ap.add_argument("--baroscale-period", type=positive_int, default=None,
                    help="DYNAMICS baroscale_period override for NPT work inputs (default: use input value)")
    ap.add_argument("--timeout", type=int, default=7200, help="per-run timeout seconds (default 7200)")
    ap.add_argument("--lock", default="/tmp/bench.lock", help="advisory serialisation lockfile")
    ap.add_argument("--allow-failures", action="store_true",
                    help="write a partial CSV and exit 0 even if selected cells fail")
    ap.add_argument("--timestamp", default=None, help="results filename stamp (else auto)")
    ap.add_argument("--out", default=None, help="explicit CSV output path (overrides --timestamp)")
    ap.add_argument("--spdyn", default=SPDYN,
                    help="path to the spdyn binary (default: ../src/spdyn_singlempi/spdyn)")
    args = ap.parse_args()
    validate_eneout_period(args, ap)

    # Use the requested spdyn binary everywhere (run_spdyn reads the module global).
    SPDYN = os.path.abspath(os.path.expanduser(args.spdyn))
    if not (os.path.isfile(SPDYN) and os.access(SPDYN, os.X_OK)):
        sys.exit("error: spdyn not found or not executable: %s" % SPDYN)

    # Resolve tune set.
    if args.full_autotune:
        tune = {"kernel", "cell", "nblist"}
    else:
        raw = [t.strip().lower() for t in args.tune.split(",") if t.strip()]
        if "all" in raw:
            tune = {"kernel", "cell", "nblist"}
        elif "none" in raw or not raw:
            tune = set()
        else:
            alias = {"nbupdate": "nblist", "nb": "nblist", "cell_size": "cell"}
            tune = set(alias.get(t, t) for t in raw)
            bad = tune - {"kernel", "cell", "nblist"}
            if bad:
                ap.error("unknown tuner(s): %s (use kernel,cell,nblist,all,none)" % ",".join(bad))
    args.tune_set = tune
    if "nblist" in args.tune_set and args.nbupdate_period is not None:
        ap.error("--nbupdate-period cannot be combined with nblist autotune; "
                 "pairlistdist and nbupdate_period are tuned as a coupled pair")

    systems = [s.strip() for s in args.systems.split(",") if s.strip()]
    ensembles = [e.strip().lower() for e in args.ensembles.split(",") if e.strip()]
    dts = []
    for d in args.dt.split(","):
        d = d.strip().replace("fs", "")
        if d in ("2", "4"):
            dts.append(d + "fs")
    if not dts:
        ap.error("--dt must include 2 and/or 4")
    for s in systems:
        if s not in ALL_SYSTEMS:
            ap.error("unknown system %r (choices: %s)" % (s, ",".join(ALL_SYSTEMS)))

    if not os.path.isfile(SPDYN):
        ap.error("spdyn binary not found at %s (build it: cd src/spdyn_singlempi && make)" % SPDYN)

    results = []
    failures = []
    csv_path = None
    log_root = None
    log("waiting for benchmark lock: %s" % args.lock)
    lock_fd = acquire_lock(args.lock, "run_benchmark")
    try:
        try:
            extracted = ensure_system_data(systems)
            stamp = args.timestamp or datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
            requested_csv_path = args.out or os.path.join(RESULTS, stamp + ".csv")
            csv_path, log_root = unique_run_paths(requested_csv_path)
            args.log_root = log_root
            args.input_log_dir = os.path.join(log_root, "inputs")
            os.makedirs(os.path.join(log_root, "autotune"), exist_ok=True)
            os.makedirs(os.path.join(log_root, "production"), exist_ok=True)
            os.makedirs(args.input_log_dir, exist_ok=True)
            open_benchmark_log(os.path.join(log_root, "benchmark.log"))
            input_columns = planned_input_columns(systems, ensembles, dts, args)
            create_results_csv(csv_path, input_columns)

            log("=== GENESIS GPU benchmark ===")
            log("spdyn      : %s" % SPDYN)
            log("systems    : %s" % ",".join(systems))
            log("ensembles  : %s" % ",".join(ensembles))
            log("dt         : %s" % ",".join(dts))
            log("tuners     : %s" % (",".join(sorted(tune)) or "none"))
            log("warmup/measure: %d/%d   nsteps meas/tune: %d/%d" %
                (args.warmup, args.measure, args.nsteps, args.tune_nsteps))
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
            for sysname in extracted:
                log("data       : extracted data/%s.tgz -> data/%s/" % (sysname, sysname))
            log("")

            for s in systems:
                for e in ensembles:
                    for d in dts:
                        r = run_cell(s, e, d, args)
                        if r:
                            results.append(r)
                            append_result_csv(csv_path, r, tune, input_columns)
                        else:
                            failures.append("%s_%s_%s" % (s, e, d))
        finally:
            release_lock(lock_fd)

        # Print and persist the aggregate table.
        table_lines = summary_table_lines(results, tune, args.measure, args.nsteps)
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
