#!/usr/bin/env python3
"""Checks of the analysis tools on synthetic data with known answers.

Usage:
    python synthetic.py                # in a new temporary folder (printed at the end)
    python synthetic.py --dir checks

Run it after installing the tools in a new environment or after an OVITO
update. It writes the inputs, runs every tool as from the command line and
prints expected vs obtained, PASS or FAIL. The OVITO tools are skipped if
ovito is not installed. Takes under a minute.

Cases:
  friction   three sliding logs, tau = 0.1 P (P = 1, 2, 5 GPa) plus
             correlated noise of known integrated correlation time (19 rows;
             the estimate scatters by about 35 % with ~170 independent
             samples); a log in units real
  restarts   a run continued from a restart written before the end of the
             first job, then a reset_timestep 0 (logs and dumps)
  vprofile   tanh profile (z0 = 20 A, w = 1.5 A: w10-90 = 2.20 w = 3.30 A)
             and linear profile between z = 10 and 30 A (Couette: rate
             0.25 x 10^10 s^-1, w10-90 = 0.80 h = 16 A)
  zprofile   atoms uniform in 10 < z < 30 A: number density N / (lx ly 20);
             an element present in half of the frames
  structure  diamond: 4 neighbours, C-C-C 109.47 deg
  bonds      C / Si+N sides with O bridges, counted by hand
  detached   a diamond block that loses an atom pair and an atom, then
             gets the atom back; an atom held to the block only through
             an O, with and without --connect
  species    two H2O, an OH bonded to an Fe atom, a lone H; an element that
             appears only in the last frame
"""
import argparse
import importlib.util
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

TOOLS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools")
results = []


def run(tool, *args):
    cmd = [sys.executable, os.path.join(TOOLS, tool), *map(str, args)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        print(p.stdout[-2000:], p.stderr[-2000:])
        raise SystemExit(f"{tool} failed")
    return p.stdout


def check(case, what, expected, got, tol):
    ok = abs(got - expected) <= tol
    results.append(ok)
    print(f"{case:>10}  {what:<28} expected {expected:>9.4g}  got {got:>9.4g}  "
          f"(tol {tol:.2g})  {'PASS' if ok else 'FAIL'}")


def write_dump(fn, frames, box, elements=None):
    """frames: list of (step, positions); box: [(lo, hi)] * 3."""
    with open(fn, "w") as f:
        for step, pos in frames:
            f.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(pos)}\n"
                    "ITEM: BOX BOUNDS pp pp pp\n" + "".join(f"{lo} {hi}\n" for lo, hi in box)
                    + "ITEM: ATOMS id element xu yu zu\n")
            el = elements if elements is not None else ["C"] * len(pos)
            f.write("".join(f"{k + 1} {e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n"
                            for k, (e, p) in enumerate(zip(el, pos))))


# ---------------------------------------------------------------- friction
def friction(rng):
    area, phi, n = 400.0, 0.9, 4000                 # tau_int = (1 + phi) / (1 - phi) = 19
    logs = []
    for P in (1, 2, 5):
        os.makedirs(f"P{P}", exist_ok=True)
        x = np.zeros(n)
        e = rng.normal(size=n)
        for k in range(1, n):
            x[k] = phi * x[k - 1] + np.sqrt(1 - phi ** 2) * e[k]
        tau = 0.1 * P + 0.03 * x                     # GPa
        fx = -tau * area / 160.21766208             # friction = -f_load[1]
        fz = P * area / 160.21766208 + 0.02 * rng.normal(size=n)
        with open(f"P{P}/slide.log", "w") as f:
            f.write("units metal\ntimestep 0.001\nrun 400000\n   Step Time f_load[1] f_load[3]\n")
            f.write("".join(f"{100 * k} {0.1 * k:.3f} {fx[k]:.8f} {fz[k]:.8f}\n" for k in range(n)))
            f.write("Loop time of 1 on 1 procs for 400000 steps with 100 atoms\n")
        logs.append(f"P{P}/slide.log")
    run("friction.py", *logs, "--force", "-f_load[1]", "--normal", "f_load[3]", "--area", area,
        "--labels", "1GPa", "2GPa", "5GPa", "-o", "friction")
    d = pd.read_csv("friction.csv")
    for (_, r), P in zip(d.iterrows(), (1, 2, 5)):
        check("friction", f"tau at {P} GPa (GPa)", 0.1 * P, r["tau_GPa"], 3 * r["tau_sem"])
    check("friction", "mu at 5 GPa", 0.1, d["mu"].iloc[-1], 3 * d["mu_err"].iloc[-1])
    check("friction", "tau_int (rows, mean of 3)", 19, d["tau_int_rows"].mean(), 8)
    # units real: forces in kcal/mol/A, tau = F/A x 6.9477 GPa
    with open("real.log", "w") as f:
        f.write("units real\ntimestep 1.0\nrun 4000\n   Step Time f_load[1]\n")
        f.write("".join(f"{k} {k:.1f} {-0.5 * area / 6.9477 * (1 + 0.01 * np.sin(k)):.8f}\n"
                        for k in range(n)))
        f.write("Loop time of 1 on 1 procs for 4000 steps with 100 atoms\n")
    run("friction.py", "real.log", "--force", "-f_load[1]", "--area", area, "-o", "friction_real")
    check("friction", "tau, units real (GPa)", 0.5, pd.read_csv("friction_real.csv")["tau_GPa"][0], 0.001)


def restarts():
    """A run continued from a restart written before the end of the first
    job (steps 1000-1500 written twice), then a reset_timestep 0."""
    def log(fn, steps):
        with open(fn, "w") as f:
            f.write("units metal\ntimestep 0.001\nrun 1\n   Step Time Temp\n")
            f.write("".join(f"{s} {s * 0.001:.3f} {300 + s / 100:.2f}\n" for s in steps))
            f.write("Loop time of 1 on 1 procs\n")
    log("a.log", range(0, 1501, 100))
    log("b.log", range(1000, 2001, 100))
    log("c.log", range(0, 501, 100))
    run("thermo.py", "a.log", "b.log", "c.log", "-c", "Temp", "-o", "restarts")
    d = pd.read_csv("restarts.csv")
    check("restarts", "log rows (21 + 6)", 27, len(d), 0)
    check("restarts", "end of stage 0 (ps)", 2.0, d[d["stage"] == 0]["time"].max(), 1e-9)
    check("restarts", "end of stage 1 (ps)", 2.5, d["time"].max(), 1e-9)
    pos = np.zeros((1, 3))
    write_dump("ra.lammpstrj", [(s, pos) for s in range(0, 1501, 500)], [(0, 10)] * 3)
    write_dump("rb.lammpstrj", [(s, pos) for s in range(1000, 2001, 500)], [(0, 10)] * 3)
    sys.path.insert(0, TOOLS)
    import lmp_io as L
    idx = L.index_dump(["ra.lammpstrj", "rb.lammpstrj"])
    check("restarts", "dump frames (0-2000 by 500)", 5, len(idx), 0)


# ---------------------------------------------------------------- vprofile
def vprofile(rng):
    n, lz = 3000, 40.0
    z = rng.uniform(0, lz, n)
    x0, y = rng.uniform(0, 20, n), rng.uniform(0, 20, n)
    box = [(0, 20), (0, 20), (0, lz)]
    for name, v in (("tanh", 0.025 * (1 + np.tanh((z - 20) / 1.5))),
                    ("couette", 0.05 * np.clip((z - 10) / 20, 0, 1))):
        frames = [(100 * k, np.c_[x0 + v * 0.1 * k + rng.normal(0, 0.005, n), y, z])
                  for k in range(60)]
        write_dump(f"{name}.lammpstrj", frames, box)
    out = run("vprofile.py", "tanh.lammpstrj", "--dt", 0.001, "--bin", 0.2, "-o", "vp_tanh")
    z0, w = map(float, re.search(r"z0 = ([\d.]+) .* w = ([\d.]+)", out).groups())
    check("vprofile", "tanh z0 (A)", 20.0, z0, 0.1)
    check("vprofile", "tanh w (A)", 1.5, w, 0.1)
    check("vprofile", "tanh w10-90 (A)", 3.30, float(re.search(r"w10-90 = ([\d.]+)", out)[1]), 0.2)
    out = run("vprofile.py", "couette.lammpstrj", "--dt", 0.001, "--bin", 0.5, "-o", "vp_couette")
    rate = float(re.search(r"shear rate ([\d.]+)", out)[1])
    check("vprofile", "Couette rate (10^10 s^-1)", 0.25, rate, 0.01)
    check("vprofile", "Couette h (A)", 20.0, float(re.search(r"h = ([\d.]+)", out)[1]), 1.0)
    check("vprofile", "Couette w10-90 (A)", 16.0, float(re.search(r"w10-90 = ([\d.]+)", out)[1]), 0.5)


# ---------------------------------------------------------------- zprofile
def zprofile(rng):
    n = 8000
    frames = [(10 * k, np.c_[rng.uniform(0, 20, n), rng.uniform(0, 20, n), rng.uniform(10, 30, n)])
              for k in range(10)]
    write_dump("slab.lammpstrj", frames, [(0, 20), (0, 20), (0, 40)])
    run("zprofile.py", "slab.lammpstrj", "--bin", 1.0, "-o", "zp")
    # an element present in half of the frames counts half
    el = ["C"] * 4000 + ["H"] * 4000
    write_dump("half.lammpstrj", [(k, frames[k][1]) if k % 2 else (k, frames[k][1][:4000])
                                  for k in range(4)], [(0, 20), (0, 20), (0, 40)], el)
    run("zprofile.py", "half.lammpstrj", "--bin", 1.0, "-o", "zp_half")
    h = pd.read_csv("zp_half.csv")
    check("zprofile", "H in 2 of 4 frames (A^-3)", 0.5 * 4000 / 8000,
          h[(h["z"] > 12) & (h["z"] < 28)]["H (A^-3)"].mean(), 0.02)
    d = pd.read_csv("zp.csv")
    mid = d[(d["z"] > 12) & (d["z"] < 28)].iloc[:, 1]
    check("zprofile", "density 12-28 A (A^-3)", n / (20 * 20 * 20), mid.mean(), 0.02 * n / 8000)
    check("zprofile", "density outside (A^-3)", 0.0, d[(d["z"] < 9) | (d["z"] > 31)].iloc[:, 1].max(), 1e-12)


# ---------------------------------------------------------------- OVITO tools
def diamond():
    from ase.build import bulk
    a = bulk("C", "diamond", a=3.567, cubic=True).repeat((3, 3, 3))
    return a.positions, a.cell.lengths()


def structure(rng):
    p, L = diamond()
    frames = [(0, p), (100, p + rng.normal(0, 0.02, p.shape))]
    write_dump("diamond.lammpstrj", frames, [(0, x) for x in L])
    out = run("structure.py", "diamond.lammpstrj", "--last", 2, "--angles", "C-C-C", "-o", "dia")
    check("structure", "mean C neighbours", 4.0, float(re.search(r"mean ([\d.]+)", out)[1]), 0.001)
    check("structure", "C-C-C mean angle (deg)", 109.47,
          float(re.search(r"C-C-C: .* mean ([\d.]+)", out)[1]), 0.3)


def bonds():
    atoms = [("C", 0, 0, 0), ("Si", 0, 0, 1.9),                         # C-Si
             ("C", 5, 0, 0), ("O", 5, 0, 1.4), ("Si", 5, 0, 3.0),       # C-O-Si
             ("C", 10, 0, 0), ("C", 11.4, 0, 1.4), ("O", 10, 0, 1.4),   # O with two C: 2 x C-O-Si
             ("Si", 10, 0, 3.0),                                        # and C(11.4)-Si at 2.13 A
             ("C", 15, 0, 0), ("N", 15, 0, 1.45), ("H", 15, 5, 0)]      # C-N
    el = [a[0] for a in atoms]
    p = np.array([a[1:] for a in atoms]) + 2.0
    write_dump("bridge.lammpstrj", [(0, p), (10, p)], [(0, 20)] * 3, el)
    run("bonds.py", "bridge.lammpstrj", "--side1", "C", "--side2", "Si+N", "--bridge", "O", "-o", "bridge")
    d = pd.read_csv("bridge.csv").iloc[-1]
    for col, exp in (("C-Si", 2), ("C-N", 1), ("C-O-Si", 3), ("bridging O", 2), ("total", 6)):
        check("bonds", col, exp, d[col], 0)
    check("bonds", "total per nm^2", 6 / 4.0, d["total_nm2"], 1e-9)


def detached():
    p0, L = diamond()
    nb = int(np.argsort(np.linalg.norm(p0 - p0[0], axis=1))[1])
    p1 = p0.copy()
    p1[[0, nb]] += [0, 0, 15.0]                 # a bonded pair leaves
    p1[150] += [0, 0, 18.0]                     # and a single atom
    p2 = p1.copy()
    p2[150] -= [0, 0, 18.0]                     # the single atom comes back
    write_dump("wear.lammpstrj", [(0, p0), (100, p1), (200, p2)],
               [(0, L[0]), (0, L[1]), (-2, 40)])
    run("detached.py", "wear.lammpstrj", "--track", "C", "-o", "wear")
    d = pd.read_csv("wear.csv")
    check("detached", "detached C, frame 1", 3, d["detached C"].iloc[1], 0)
    check("detached", "fragments, frame 1", 2, d["fragments"].iloc[1], 0)
    check("detached", "detached C, frame 2", 2, d["detached C"].iloc[2], 0)
    f = pd.read_csv("wear_fragments.csv")
    ok = list(f["formula"]) == ["C2"]
    results.append(ok)
    print(f"{'detached':>10}  {'fragment formula':<28} expected        C2  got {','.join(f['formula']):>9}"
          f"  {'PASS' if ok else 'FAIL'}")
    # a C atom held to the block only through an O (O counts with --connect)
    top = int(np.argmax(p0[:, 2]))
    o = p0[top] + [0, 0, 1.4]
    far = p0[top] + [0, 0, 15]
    frames = [(0, np.vstack([p0, p0[top] + [0, 0, 1.5], far])),
              (100, np.vstack([p0, o + [1.21, 0, 0.7], o]))]
    write_dump("held.lammpstrj", frames, [(0, L[0]), (0, L[1]), (-2, 40)], ["C"] * len(p0) + ["C", "O"])
    run("detached.py", "held.lammpstrj", "--track", "C", "-o", "held")
    check("detached", "held by O, followed only", 1, pd.read_csv("held.csv")["detached"].iloc[1], 0)
    run("detached.py", "held.lammpstrj", "--track", "C", "--connect", "C+O", "-o", "held_c")
    check("detached", "held by O, --connect C+O", 0, pd.read_csv("held_c.csv")["detached"].iloc[1], 0)


def species():
    def h2o(o):
        o = np.array(o, float)
        return [("O", *o), ("H", *(o + [0.757, 0.586, 0])), ("H", *(o + [-0.757, 0.586, 0]))]
    atoms = h2o([3, 3, 3]) + h2o([10, 3, 3]) + [
        ("Fe", 3, 10, 3), ("O", 3, 10, 5.0), ("H", 3, 10, 5.97),         # OH on Fe (Fe-O 2.0 A)
        ("H", 12, 12, 12)]
    el = [a[0] for a in atoms]
    p = np.array([a[1:] for a in atoms], float)
    write_dump("water.lammpstrj", [(0, p), (10, p)], [(0, 20)] * 3, el)
    run("species.py", "water.lammpstrj", "--track", "H+O", "-o", "water")
    # H appears only in the last frame (deposition): no crash, counted
    q = np.array([[3, 10, 3], [3, 10, 5.0], [3, 10, 5.97], [12, 12, 12]])
    write_dump("late.lammpstrj", [(0, q[:2]), (10, q)], [(0, 20)] * 3, ["Fe", "O", "H", "H"])
    run("species.py", "late.lammpstrj", "--track", "H+O", "-o", "late")
    check("species", "lone H, element added later", 1, pd.read_csv("late.csv").iloc[-1].get("H", -1), 0)
    d = pd.read_csv("water.csv").iloc[-1]
    for col, exp in (("H2O", 2), ("HO (bound)", 1), ("H", 1), ("O bonded to untracked", 1)):
        check("species", col, exp, d.get(col, -1), 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", help="working folder (default: a new temporary folder)")
    a = ap.parse_args()
    if a.dir:
        os.makedirs(a.dir, exist_ok=True)
    os.chdir(a.dir or tempfile.mkdtemp(prefix="synthetic_check_"))
    rng = np.random.default_rng(7)
    friction(rng)
    restarts()
    vprofile(rng)
    zprofile(rng)
    if importlib.util.find_spec("ovito") is None:
        print("ovito not installed: structure, bonds, detached, species skipped")
    else:
        structure(rng)
        bonds()
        detached()
        species()
    print(f"\n{sum(results)}/{len(results)} checks passed; files in {os.getcwd()}")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
