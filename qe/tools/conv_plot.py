#!/usr/bin/env python3
"""Convergence table and plot from a series of pw.x outputs.

Usage:
    python conv_plot.py ecut_*.pwo            # x value read from the file name
    python conv_plot.py k_*.pwo --tol 1.0     # tolerance in meV/atom
    python conv_plot.py degauss_*.pwo --ref min

The reference is the converged end of the scan: the largest value for cutoffs
and k-grids (default, --ref max), the smallest for the smearing width, which
converges for degauss -> 0 (--ref min). For each run it prints the energy per
atom, the difference from the reference in meV/atom, the pressure (if tstress)
and the total force, then saves conv.png. It reports the range of values whose
|dE| stays below --tol.
Needs numpy and matplotlib.
"""
import argparse
import re
import sys

import numpy as np

from pwo_summary import PATTERNS, RY_TO_EV, last


def read_point(path):
    text = open(path, errors="replace").read()
    if "JOB DONE." not in text or "convergence NOT achieved" in text:
        print(f"warning: {path} not finished or SCF not converged, skipped")
        return None
    energy = float(last(PATTERNS["energy"], text))
    nat = int(last(PATTERNS["nat"], text))
    press = last(PATTERNS["press"], text)
    force = last(PATTERNS["force"], text)
    match = re.search(r"([A-Za-z]+)_(\d+(?:\.\d+)?)\.pwo$", path)
    if not match:
        sys.exit(f"cannot read the parameter value from the name {path}")
    return {
        "name": match.group(1),
        "x": float(match.group(2)),
        "e": energy * RY_TO_EV / nat,
        "p": float(press) if press else np.nan,
        "f": float(force) if force else np.nan,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--tol", type=float, default=1.0, help="meV/atom")
    parser.add_argument("--out", default="conv.png")
    parser.add_argument("--ref", choices=("max", "min"), default="max",
                        help="reference = largest (cutoff, k) or smallest (degauss) value")
    args = parser.parse_args()

    points = [p for p in (read_point(f) for f in args.files) if p]
    if len(points) < 2:
        sys.exit("need at least two finished runs")
    points.sort(key=lambda p: p["x"])
    x = np.array([p["x"] for p in points])
    ref = points[-1] if args.ref == "max" else points[0]
    de = (np.array([p["e"] for p in points]) - ref["e"]) * 1000.0
    dp = np.array([p["p"] for p in points]) - ref["p"]

    print(f"{'x':>10s} {'E/atom (eV)':>14s} {'dE (meV/at)':>12s} {'P (kbar)':>10s} "
          f"{'dP (kbar)':>10s} {'F_tot (Ry/b)':>12s}")
    for p, d, dpp in zip(points, de, dp):
        print(f"{p['x']:10g} {p['e']:14.6f} {d:12.3f} {p['p']:10.2f} {dpp:10.2f} {p['f']:12.5f}")

    ok = np.abs(de) < args.tol
    if args.ref == "max":
        first = next((x[i] for i in range(len(x)) if ok[i:].all()), None)
        print(f"\n|dE| < {args.tol} meV/atom from x = {first:g} onwards (reference = largest x)")
    else:
        lastok = next((x[i] for i in range(len(x) - 1, -1, -1) if ok[:i + 1].all()), None)
        print(f"\n|dE| < {args.tol} meV/atom up to x = {lastok:g} (reference = smallest x)")

    labels = {"ecut": "ecutwfc (Ry)", "ecutrho": "ecutrho (Ry)", "k": "k-points per direction",
              "degauss": "degauss (Ry)", "scale": "lattice scale factor"}
    xlabel = labels.get(points[0]["name"], points[0]["name"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_p = not np.all(np.isnan(dp))
    fig, axes = plt.subplots(1, 2 if has_p else 1, figsize=(9 if has_p else 4.5, 3.5),
                             squeeze=False)
    ax = axes[0, 0]
    ax.plot(x, de, "o-")
    ax.axhspan(-args.tol, args.tol, color="0.85")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("E - E_ref (meV/atom)")
    ax.axvline(ref["x"], color="0.5", lw=0.8, ls=":")
    if has_p:
        ax = axes[0, 1]
        ax.plot(x, dp, "o-")
        ax.axhspan(-0.5, 0.5, color="0.85")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("P - P_ref (kbar)")
    for a in axes[0]:
        a.xaxis.set_major_locator(plt.MaxNLocator(6))
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("saved", args.out)


if __name__ == "__main__":
    main()
