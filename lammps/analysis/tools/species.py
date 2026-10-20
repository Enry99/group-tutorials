#!/usr/bin/env python3
"""Molecule census vs time (OVITO bonds, clusters with scipy): the
molecules and fragments formed by a set of atoms (lubricant, water,
additive), free or bound to the rest of the system (e.g. a surface).

Usage:
    python species.py highT.lammpstrj --track molecules --groups Fe_water.groups.lmp \\
        --cutoff Fe-Fe:2.7 --dt 0.0005
    python species.py dyn*.lammpstrj --track C+H+O --dt 0.0005 --every 10 --top 8

What it does:
  1. tracked atoms: elements, ids or groups of a groups file, joined
     with '+' (as detached.py);
  2. bonds with the pairwise cutoffs of cutoffs.py (--cutoff A-B:R,
     --alternative; check them with structure.py: the free/bound split
     depends on the cutoff between tracked and untracked atoms, e.g.
     Fe-O); molecules = clusters of tracked atoms connected by bonds;
     a molecule with a bond to an untracked atom is 'bound';
  3. per frame: number of molecules of each formula (Hill order), free
     and bound, and the tracked atoms bonded to untracked ones, per
     element;
  4. frames: all by default (--every, --last, --tmin/--tmax with --dt;
     restart duplicates dropped); writes <out>.csv (one column per
     species) and <out>.png (the --top species by mean count; --rolling
     N as in bonds.py).
Needs ovito, numpy, pandas, matplotlib, scipy, ase.
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ovito_io as O  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#c43c6a", "#5a5a5a", "#9c6b30"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    O.add_common_args(parser)
    parser.add_argument("--track", required=True, help="elements 'C+H+O', ids '1:500', groups 'g1+g2'")
    parser.add_argument("--groups", help="groups file for group names")
    parser.add_argument("--top", type=int, default=6, help="species in the plot (default 6)")
    parser.add_argument("--rolling", type=int, help="frames of the centred rolling mean")
    parser.add_argument("-o", "--output", default="species")
    args = parser.parse_args()

    pipeline, steps, els = O.load(args.dumps, args.type_names)
    track, desc = O.selector(args.track, els, args.groups)
    frames, t, stage = O.window(steps, args)
    O.add_bonds(pipeline, els, args)
    print(f"tracked: {desc}\n{len(frames)} frames, steps {steps[frames[0]]}-{steps[frames[-1]]}")

    rows = []
    for k, tk, sk in zip(frames, t, stage):
        data = pipeline.compute(k)
        el, ids = O.elements_of(data), O.ids_of_frame(data)
        i, j, _ = O.bond_list(data)
        f = track(el, ids)
        lab = O.clusters(len(el), i, j, f)
        out = f[i] != f[j]                                  # tracked - untracked bonds
        anchored = np.zeros(len(el), dtype=bool)
        anchored[np.where(f[i[out]], i[out], j[out])] = True
        bound = set(lab[anchored])
        census = Counter(fm + (" (bound)" if c in bound else "")
                         for c, fm in O.formulas(lab, el).items())
        row = {"frame": k, "step": steps[k], "time": tk, "stage": sk, "molecules": int(lab.max() + 1)}
        row.update(census)
        if not rows and not f.any():
            sys.exit("no tracked atom in the first frame: check --track")
        for e in sorted(set(el[f]), key=lambda e: (e not in els, els.index(e) if e in els else 0)):
            row[f"{e} bonded to untracked"] = int(np.count_nonzero(anchored & (el == e)))
        rows.append(row)
    df = pd.DataFrame(rows).fillna(0)
    meta = ["frame", "step", "time", "stage", "molecules"]
    anch = [c for c in df.columns if c.endswith("bonded to untracked")]
    spec = sorted([c for c in df.columns if c not in meta + anch], key=lambda c: -df[c].mean())
    df[spec + anch] = df[spec + anch].astype(int)
    df = df[meta + spec + anch]
    df.to_csv(args.output + ".csv", index=False)

    print(f"\n{'species':>20} {'first':>7} {'last':>7} {'mean':>8}")
    for c in spec[:20]:
        print(f"{c:>20} {df[c].iloc[0]:>7} {df[c].iloc[-1]:>7} {df[c].mean():>8.1f}")
    if len(spec) > 20:
        print(f"{'':>20} ... {len(spec) - 20} more in the csv")
    for c in anch:
        print(f"{c:>20}: first {df[c].iloc[0]}, last {df[c].iloc[-1]}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for k, c in enumerate(spec[:min(args.top, len(COLORS))]):
        O.plot_series(ax, df["time"], df[c], COLORS[k], c, args.rolling)
    O.stage_lines(ax, df["time"], df["stage"])
    ax.set(xlabel=O.time_label(args), ylabel="molecules")
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    print(f"-> {args.output}.csv, {args.output}.png")


if __name__ == "__main__":
    main()
