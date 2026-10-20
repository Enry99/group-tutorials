"""Helpers shared by the analysis scripts (not meant to be run): LAMMPS
logs (thermo tables), text dumps (frames sorted by id), groups files of
structures/tools/interface.py, time axis across reset_timestep, and the
statistics of correlated series (blocking, integrated correlation time).

Readers without OVITO: numpy and pandas only. The blocking and tau_int
functions are those of friction_vs_load.py (E. Pedretti).
"""
import glob
import io
import mmap
import re
import sys

import numpy as np
import pandas as pd

EV_PER_A3_TO_GPA = 160.21766208         # 1 eV/A^3 in GPa
TIME_UNIT = {"metal": "ps", "real": "fs"}


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def expand(patterns):
    """Glob patterns and/or file names -> list in natural order (1, 2, 10),
    patterns expanded one by one, so the order of the arguments is kept."""
    if isinstance(patterns, str):
        patterns = [patterns]
    out = []
    for p in patterns:
        hits = sorted(glob.glob(p), key=natural_key)
        if not hits:
            raise SystemExit(f"no file matches '{p}'")
        out += hits
    return out


# ---------------------------------------------------------------- logs
def read_log(path):
    """Thermo tables of a LAMMPS log, one dict per run:
    {"columns": {name: array}, "dt": timestep, "units": str, "log": path,
    "minimize": bool}. WARNING lines inside a table are skipped; the table
    ends at "Loop time" or at the first line that is not a row (e.g. an
    ERROR). A table after a 'minimize' command is marked as minimization."""
    runs, dt, units = [], None, None             # units: None if the log has no 'units' (restart)
    minimize = False
    lines = open(path, errors="replace").read().splitlines()
    i = 0
    while i < len(lines):
        w = lines[i].split()
        if w[:1] == ["minimize"]:
            minimize = True
        elif w[:1] == ["run"]:
            minimize = False
        if w[:1] == ["units"] and len(w) > 1:
            units = w[1]
        elif w[:1] == ["timestep"] and len(w) > 1:
            try:
                dt = float(w[1])
            except ValueError:               # unsubstituted variable: next echo
                pass
        elif w[:1] == ["Step"] or w[:2] == ["Step", "Time"]:
            head, rows = w, []
            i += 1
            while i < len(lines) and not lines[i].startswith("Loop time"):
                v = lines[i].split()
                if lines[i].startswith("WARNING") or not v:
                    i += 1
                    continue
                if len(v) != len(head):
                    break
                try:
                    rows.append([float(x) for x in v])
                except ValueError:
                    break
                i += 1
            if rows:
                a = np.array(rows)
                runs.append({"columns": {h: a[:, k] for k, h in enumerate(head)},
                             "dt": dt, "units": units, "log": path, "minimize": minimize})
            continue
        i += 1
    return runs


def join_runs(runs, dts=None, minimize=False):
    """Concatenate thermo tables (runs of one or several logs, in order) on
    one time axis. A new stage starts where the step counter goes back
    (reset_timestep) or a new log begins; stages are laid end to end. The
    first row of a run that repeats the last row of the previous one is
    dropped. Time: the Time column if printed, else Step x timestep, run by
    run (a change of timestep inside a stage is followed).
    dts: timestep per stage (time units of the log), to override the logs.
    Minimizations are left out unless minimize=True (their steps are
    iterations, not time). Returns a DataFrame with the columns, 'time'
    and 'stage'."""
    out, stage, t_end, prev = [], -1, 0.0, None
    for r in runs:
        if r["minimize"] and not minimize:
            continue
        d = pd.DataFrame(r["columns"])
        step = d["Step"].to_numpy()
        new_stage = prev is None or step[0] < prev["step"] or r["log"] != prev["log"]
        if new_stage:
            stage += 1
        dt = dts[min(stage, len(dts) - 1)] if dts else r["dt"]
        if "Time" in d:
            tc = d["Time"].to_numpy()
            ref = tc[0] if new_stage else prev["time_col"]
        elif dt is not None:
            tc = step * dt
            ref = tc[0] if new_stage else prev["step"] * dt
        else:
            raise SystemExit(f"{r['log']}: no Time column and no timestep: give --dt")
        d = d.assign(time=t_end + tc - ref, stage=stage)
        if not new_stage and step[0] == prev["step"]:
            d = d.iloc[1:]
        if d.empty:
            continue
        out.append(d)
        t_end = d["time"].iloc[-1]
        prev = {"step": step[-1], "log": r["log"], "time_col": tc[-1]}
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------- dumps
def index_dump(files, dedup=True):
    """(file, byte offset, step) of every frame of text dumps, in file order.
    A frame with the same step as the frame before it (restart overlap: the
    last frame of a file written again at the start of the next) is dropped,
    the later copy kept (dedup=False keeps both). Steps that come back after
    a reset_timestep are different frames and are all kept."""
    out = []
    for fn in files:
        with open(fn, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            pos = mm.find(b"ITEM: TIMESTEP")
            while pos >= 0:
                nl = mm.find(b"\n", pos)
                step = int(mm[nl + 1:mm.find(b"\n", nl + 1)])
                if dedup and out and out[-1][2] == step:
                    out.pop()
                out.append((fn, pos, step))
                pos = mm.find(b"ITEM: TIMESTEP", pos + 14)
    return out


def read_frame(fn, pos):
    """One frame: {"step", "lo", "hi", "tilt", "pbc", "cols": DataFrame
    sorted by id if present}."""
    with open(fn, "rb") as fh:
        fh.seek(pos)
        fh.readline()
        step = int(fh.readline())
        fh.readline()
        n = int(fh.readline())
        hdr = fh.readline().decode().split()               # ITEM: BOX BOUNDS [xy xz yz] pp pp pp
        b = [list(map(float, fh.readline().split())) for _ in range(3)]
        cols = fh.readline().decode().split()[2:]
        block = b"".join(fh.readline() for _ in range(n))
    tilt = [r[2] if len(r) > 2 else 0.0 for r in b]
    lo = np.array([r[0] for r in b])
    hi = np.array([r[1] for r in b])
    if "xy" in hdr:                                       # bounding box -> box (LAMMPS manual, dump)
        xy, xz, yz = tilt
        lo[0] -= min(0.0, xy, xz, xy + xz)
        hi[0] -= max(0.0, xy, xz, xy + xz)
        lo[1] -= min(0.0, yz)
        hi[1] -= max(0.0, yz)
    d = pd.read_csv(io.BytesIO(block), sep=r"\s+", header=None, names=cols)
    if "id" in d:
        d = d.sort_values("id", kind="stable").reset_index(drop=True)
    return {"step": step, "lo": lo, "hi": hi, "tilt": np.array(tilt),
            "pbc": [x == "pp" for x in hdr[-3:]], "cols": d}


def frames(files, first=0, last=None, every=1):
    """Iterate over the frames of text dumps (list or patterns), with the
    frame index: yields (k, frame)."""
    idx = index_dump(expand(files))
    sel = list(range(len(idx)))[first:last][::every]
    for k in sel:
        fn, pos, _ = idx[k]
        yield k, read_frame(fn, pos)


def window(steps, dt=None, tmin=None, tmax=None, last=None, every=1):
    """Frames to analyse: those with tmin <= t <= tmax (time from
    stage_time, needs dt), then the last `last` of them, then every
    `every`-th ending at the last one. Returns (indices, t); t is the
    step count laid end to end if dt is None."""
    steps = np.asarray(steps)
    if (tmin is not None or tmax is not None) and not dt:
        raise SystemExit("--tmin/--tmax need --dt (timestep, one per stage if they differ)")
    t, _ = stage_time(steps, dt or [1.0])
    sel = np.arange(len(steps))
    sel = sel[(t >= (-np.inf if tmin is None else tmin)) & (t <= (np.inf if tmax is None else tmax))]
    if last:
        sel = sel[-last:]
    sel = sel[::-1][::every][::-1]
    if len(sel) == 0:
        raise SystemExit("no frame in the window")
    return sel, t


def stage_time(steps, dts):
    """Time of each frame from its step counter: a new stage starts where
    the counter goes back (reset_timestep), stages are laid end to end.
    dts: timestep per stage (one value: all stages). Returns (time, stage)."""
    steps = np.asarray(steps)
    stage = np.concatenate(([0], np.cumsum(steps[1:] < steps[:-1])))
    t, offset = np.zeros(len(steps)), 0.0
    for s in range(stage.max() + 1):
        m = stage == s
        dt = dts[min(s, len(dts) - 1)]
        t[m] = offset + (steps[m] - steps[m][0]) * dt
        offset = t[m][-1]
    return t, stage


# ---------------------------------------------------------------- groups
def read_groups(path):
    """'group NAME id <> lo hi' and 'group NAME union A B ...' lines of a
    groups file (interface.py) -> {name: sorted array of ids}."""
    groups = {}
    for line in open(path):
        w = line.split("#")[0].split()
        if len(w) >= 6 and w[0] == "group" and w[2] == "id" and w[3] == "<>":
            groups[w[1]] = np.arange(int(w[4]), int(w[5]) + 1)
        elif len(w) >= 4 and w[0] == "group" and w[2] == "union":
            groups[w[1]] = np.unique(np.concatenate([groups[g] for g in w[3:]]))
    return groups


def ids_of(spec, groups_file=None):
    """Atom ids from 'lo:hi', 'a,b,c' or group names of a groups file
    (several joined with '+', e.g. surf_top+thermo_top)."""
    if re.fullmatch(r"\d+:\d+", spec):
        lo, hi = map(int, spec.split(":"))
        return np.arange(lo, hi + 1)
    if re.fullmatch(r"\d+(,\d+)*", spec):
        return np.array(sorted(map(int, spec.split(","))))
    if groups_file is None:
        raise SystemExit(f"'{spec}': group names need --groups FILE")
    g = read_groups(groups_file)
    missing = [s for s in spec.split("+") if s not in g]
    if missing:
        raise SystemExit(f"groups {missing} not in {groups_file} ({', '.join(g)})")
    return np.unique(np.concatenate([g[s] for s in spec.split("+")]))


# ---------------------------------------------------------------- statistics
def blocking_sem(x, min_blocks=8, return_curve=False):
    """Flyvbjerg-Petersen blocking estimate of the standard error of the
    mean: the series is averaged in pairs again and again; the naive SEM
    grows with the block size until blocks are uncorrelated, then stays
    flat. The maximum over the levels (at least min_blocks blocks) reads
    that plateau, slightly on the high side."""
    x = np.asarray(x, dtype=float)
    sems, sizes, size = [], [], 1
    while len(x) >= min_blocks:
        sems.append(x.std(ddof=1) / np.sqrt(len(x)))
        sizes.append(size)
        if len(x) % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
        size *= 2
    sems, sizes = np.array(sems), np.array(sizes)
    return (sems.max(), sizes, sems) if return_curve else sems.max()


def tau_int(x, c=5.0):
    """Integrated autocorrelation time (in samples), Sokal window.
    N_eff = N / tau_int independent samples."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    x = x - x.mean()
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, nfft)
    acf = np.fft.irfft(f * np.conj(f), nfft)[:n].real
    if acf[0] == 0:
        return 1.0
    acf /= acf[0]
    taus = 1.0 + 2.0 * np.cumsum(acf[1:])
    w = np.arange(1, n)
    ok = w >= c * taus
    k = np.argmax(ok) if ok.any() else n - 2
    return max(float(taus[k]), 1.0)


def parse_args(parser):
    """parser.parse_args(), with expressions that start with a minus sign
    ("-f_load[1]") taken as values, not as unknown options."""
    known = set(parser._option_string_actions)
    argv = [" " + a if a.startswith("-") and a.split("=")[0] not in known
            and not re.fullmatch(r"-[\d.]+(e-?\d+)?", a) else a for a in sys.argv[1:]]
    return parser.parse_args(argv)


def evaluate(expr, df):
    """Value of an expression of thermo columns, e.g. "-f_load[1]" or
    "0.5*(f_hold_bot[1]-f_load[1])" (numpy functions as np.xxx)."""
    names = sorted(df.columns, key=len, reverse=True)
    env, code = {"np": np}, expr
    for k, c in enumerate(names):                 # whole names only (pe is not in np.exp)
        pat = r"(?<![\w.\[])" + re.escape(c) + r"(?![\w\[])"
        if re.search(pat, code):
            code = re.sub(pat, f"__c{k}", code)
            env[f"__c{k}"] = df[c].to_numpy()
    try:
        return np.asarray(eval(code.strip(), {"__builtins__": {}}, env), dtype=float)
    except Exception as e:
        raise SystemExit(f"cannot evaluate '{expr}' ({e}); columns: {', '.join(df.columns)}")


def rolling(y, window):
    """Centred rolling mean (window in samples), NaN at the edges."""
    return pd.Series(y).rolling(window, center=True).mean().to_numpy()
