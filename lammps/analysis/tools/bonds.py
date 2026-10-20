#!/usr/bin/env python3
"""Bonds across an interface vs time (OVITO): direct bonds between two
sides and bonds bridged by a third atom, per element pair, per nm^2.

Usage:
    python bonds.py 01_indent.lammpstrj --side1 sub_top --side2 sub_bot \\
        --groups Fe_tip_on_Fe.groups.lmp --cutoff Fe-Fe:2.7 --dt 0.001
    python bonds.py load.lammpstrj slide.lammpstrj --side1 C --side2 Si+N --bridge O \\
        --dt 0.0005 --every 10 --rolling 5 -o C_SiN

What it does (based on compute_bonds.py, E. Pedretti):
  1. sides: elements ('C', 'Si+N': all atoms of those elements, his
     convention when each element sits on one side only) or atoms by id
     ('1:500', groups of a groups file with --groups, joined with '+');
     atoms of the --bridge elements belong to neither side;
  2. bonds with the pairwise cutoffs of cutoffs.py (--cutoff A-B:R,
     --alternative; check them first with structure.py);
  3. per frame: direct bonds side1-side2 per element pair (A-B, A on
     side 1); bridges A-X-B, X a bridge atom bonded to A on side 1 and B
     on side 2, counted once per path (sum over X of n_A x n_B); the
     number of bridging X atoms; densities per nm^2 of interface,
     area |a x b| of the cell at that frame x --interfaces (2 for a
     periodic sandwich whose two contacts both count, 1 otherwise);
  4. frames: all by default (--every, --last, --tmin/--tmax with --dt;
     time from --dt per stage, stages laid end to end, restart duplicates
     dropped); writes <out>.csv (counts and _nm2 columns) and <out>.png
     (densities vs time, or counts with --counts, e.g. for a tip whose
     contact is smaller than the cell; stage boundaries dashed;
     --rolling N: centred mean over N frames on top of the faded data).
Needs ovito, numpy, pandas, matplotlib, scipy.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ovito_io as O  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#c43c6a"]


def count_frame(data, side1, side2, bridge_els):
    el = O.elements_of(data)
    ids = O.ids_of_frame(data)
    n = len(el)
    br = np.isin(el, bridge_els)
    s1 = side1(el, ids) & ~br
    s2 = side2(el, ids) & ~br
    if (s1 & s2).any():
        raise SystemExit(f"{np.count_nonzero(s1 & s2)} atoms are on both sides")
    i, j, _ = O.bond_list(data)
    row = {}
    d = (s1[i] & s2[j]) | (s2[i] & s1[j])
    a = np.where(s1[i[d]], el[i[d]], el[j[d]])
    b = np.where(s1[i[d]], el[j[d]], el[i[d]])
    for key, c in zip(*np.unique(np.char.add(np.char.add(a, "-"), b), return_counts=True)):
        row[str(key)] = int(c)
    for x in bridge_els:
        mx = el == x
        n_brg = np.zeros(n, dtype=bool)
        for ea in np.unique(el[s1]):
            na = n_neigh(i, j, n, mx, s1 & (el == ea))
            for eb in np.unique(el[s2]):
                nb = n_neigh(i, j, n, mx, s2 & (el == eb))
                c = int(np.sum(na * nb))
                if c:
                    row[f"{ea}-{x}-{eb}"] = c
                n_brg |= (na > 0) & (nb > 0)
        row[f"bridging {x}"] = int(n_brg.sum())
    return row, O.area(data), int(s1.sum()), int(s2.sum())


def n_neigh(i, j, n, centre, neigh):
    """For each atom in centre, the number of bonded neighbours in neigh."""
    c = np.zeros(n, dtype=np.int64)
    s = centre[i] & neigh[j]
    np.add.at(c, i[s], 1)
    s = centre[j] & neigh[i]
    np.add.at(c, j[s], 1)
    return c


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    O.add_common_args(parser)
    parser.add_argument("--side1", required=True, help="elements 'C+H', ids '1:500', groups 'g1+g2'")
    parser.add_argument("--side2", required=True)
    parser.add_argument("--groups", help="groups file for group names")
    parser.add_argument("--bridge", nargs="+", default=[], help="bridging elements, e.g. O")
    parser.add_argument("--interfaces", type=int, default=1, help="contacts in the cell (default 1)")
    parser.add_argument("--rolling", type=int, help="frames of the centred rolling mean")
    parser.add_argument("--counts", action="store_true", help="plot counts, not per nm^2 (a tip)")
    parser.add_argument("-o", "--output", default="bonds")
    args = parser.parse_args()

    pipeline, steps, els = O.load(args.dumps, args.type_names)
    side1, d1 = O.selector(args.side1, els, args.groups)
    side2, d2 = O.selector(args.side2, els, args.groups)
    for x in args.bridge:
        if x not in els:
            print(f"bridge element {x} not in the dump: no bridges")
    frames, t, stage = O.window(steps, args)
    O.add_bonds(pipeline, els, args)
    print(f"side 1: {d1}\nside 2: {d2}" + (f"\nbridges: {', '.join(args.bridge)}" if args.bridge else ""))
    print(f"{len(frames)} frames, steps {steps[frames[0]]}-{steps[frames[-1]]}")

    rows = []
    for k, tk, sk in zip(frames, t, stage):
        row, a, n1, n2 = count_frame(pipeline.compute(k), side1, side2, args.bridge)
        if not rows:
            print(f"first frame: {n1} atoms on side 1, {n2} on side 2")
            if n1 == 0 or n2 == 0:
                sys.exit("a side is empty: check --side1/--side2")
        rows.append({"frame": k, "step": steps[k], "time": tk, "stage": sk, "area_A2": a, **row})
    df = pd.DataFrame(rows).fillna(0)
    kinds = [c for c in df.columns if c not in ("frame", "step", "time", "stage", "area_A2")]
    links = [c for c in kinds if not c.startswith("bridging")]
    df[kinds] = df[kinds].astype(int)
    df["total"] = df[links].sum(axis=1) if links else 0
    per_nm2 = 100.0 / (df["area_A2"] * args.interfaces)
    for c in links + ["total"]:
        df[c + "_nm2"] = df[c] * per_nm2
    df.to_csv(args.output + ".csv", index=False)

    sfx = "" if args.counts else "_nm2"
    print(f"\n{'linkage':>12} {'first':>8} {'last':>8} {'max':>8}   "
          + ("(counts)" if args.counts else "(per nm^2; counts in the csv)"))
    for c in links + ["total"]:
        v = df[c + sfx]
        print(f"{c:>12} {v.iloc[0]:>8.4g} {v.iloc[-1]:>8.4g} {v.max():>8.4g}")
    for c in kinds:
        if c.startswith("bridging"):
            print(f"{c:>12} {df[c].iloc[0]:>8d} {df[c].iloc[-1]:>8d} {df[c].max():>8d}   (atoms)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(figsize=(6.5, 4))
    order = sorted(links, key=lambda c: -df[c].mean())
    for k, c in enumerate(order[:len(COLORS)]):
        O.plot_series(ax, df["time"], df[c + sfx], COLORS[k], c, args.rolling)
    if len(links) > 1:
        O.plot_series(ax, df["time"], df["total" + sfx], "k", "total", args.rolling, ls="--")
    if len(order) > len(COLORS):
        print(f"plot: the {len(COLORS)} most frequent kinds; all in the csv")
    O.stage_lines(ax, df["time"], df["stage"])
    ax.set(xlabel=O.time_label(args), ylabel="bonds" if args.counts else "bonds per nm$^2$")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    print(f"-> {args.output}.csv, {args.output}.png")


if __name__ == "__main__":
    main()
