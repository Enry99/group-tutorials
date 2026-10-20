#!/usr/bin/env python3
"""Atoms that detach from a body (asperity, tip, slab, coating) along a
trajectory (OVITO bonds, clusters with scipy), per element, with the
fragments they form.

Usage:
    python detached.py dyn*.lammpstrj --track C+H+O --dt 0.0005 --every 10
    python detached.py 02_slide.lammpstrj --track sub_top --groups Fe_water_Fe.groups.lmp \\
        --cutoff Fe-Fe:2.7 --dt 0.0005 -o tip_wear
    python detached.py dyn*.lammpstrj --track C+H+O --connect C+H+O --alternative   # as ovito_cho.py

What it does (based on ovito_cho.py / detached_atoms.py, E. Pedretti):
  1. tracked atoms: elements ('C+H+O'), ids ('1:500') or groups of a
     groups file (--groups), joined with '+';
  2. bonds with the pairwise cutoffs of cutoffs.py (--cutoff A-B:R,
     --alternative to reproduce results obtained with ovito_cho.py; check
     the cutoffs first with structure.py); clusters = atoms connected by
     bonds among the tracked atoms only (a bond to an untracked atom does
     not hold an atom to the body); --connect SPEC: the atoms of SPEC also
     connect clusters (--connect C+H+O as ovito_cho.py, which clusters all
     C, H, O atoms: an asperity atom held only through an adsorbed water O
     then counts as attached);
  3. reference frame = the first of the window: its clusters are printed
     (check that the largest is the body); only the atoms of the largest
     one are followed (--all-tracked: all tracked atoms);
  4. per frame: detached = followed atoms outside the cluster that holds
     most of them (the body), per element; fragments = clusters with
     detached atoms; at the last frame, the fragments by formula
     (<out>_fragments.csv; with --connect the formula includes the
     connecting atoms of the fragment);
  5. frames: all by default (--every, --last, --tmin/--tmax with --dt;
     restart duplicates dropped); writes <out>.csv and <out>.png
     (detached atoms and fragments vs time; --rolling N as in bonds.py).
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

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#c43c6a"]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    O.add_common_args(parser)
    parser.add_argument("--track", required=True, help="elements 'C+H', ids '1:500', groups 'g1+g2'")
    parser.add_argument("--groups", help="groups file for group names")
    parser.add_argument("--all-tracked", action="store_true",
                        help="follow all tracked atoms, not only the largest cluster")
    parser.add_argument("--connect", help="atoms that also connect clusters (elements, ids, groups); "
                        "'--connect C+H+O' as ovito_cho.py")
    parser.add_argument("--rolling", type=int, help="frames of the centred rolling mean")
    parser.add_argument("-o", "--output", default="detached")
    args = parser.parse_args()

    pipeline, steps, els = O.load(args.dumps, args.type_names)
    track, desc = O.selector(args.track, els, args.groups)
    connect, cdesc = O.selector(args.connect, els, args.groups) if args.connect else (None, "")
    frames, t, stage = O.window(steps, args)
    O.add_bonds(pipeline, els, args)
    print(f"tracked: {desc}" + (f"\nconnecting: {cdesc}" if args.connect else "")
          + f"\n{len(frames)} frames, steps {steps[frames[0]]}-{steps[frames[-1]]}")

    data = pipeline.compute(frames[0])
    el, ids = O.elements_of(data), O.ids_of_frame(data)
    i, j, _ = O.bond_list(data)
    m = track(el, ids)
    lab = O.clusters(len(el), i, j, m)
    form = O.formulas(lab, el)
    nclu = lab.max() + 1
    print(f"\nreference frame (step {steps[frames[0]]}): {m.sum()} tracked atoms in {nclu} cluster(s)")
    for c in range(min(nclu, 5)):
        print(f"  cluster {c + 1}: {np.count_nonzero(lab == c):>7} atoms   {form[c]}")
    if nclu > 5:
        print(f"  ... {nclu - 5} more")
    follow = ids[m] if args.all_tracked else ids[lab == 0]
    tracked_els = sorted(set(el[np.isin(ids, follow)]), key=els.index)
    print(f"-> following {len(follow)} atoms ({'all tracked' if args.all_tracked else 'cluster 1'}): "
          "check that cluster 1 is the body\n")

    rows = []
    for k, tk, sk in zip(frames, t, stage):
        data = pipeline.compute(k)
        el, ids = O.elements_of(data), O.ids_of_frame(data)
        i, j, _ = O.bond_list(data)
        f = np.isin(ids, follow)
        lab = O.clusters(len(el), i, j, f | connect(el, ids) if args.connect else f)
        body = np.bincount(lab[f]).argmax() if f.any() else -1
        det = f & (lab != body)
        row = {"frame": k, "step": steps[k], "time": tk, "stage": sk, "followed": int(f.sum())}
        for e in tracked_els:
            row[f"detached {e}"] = int(np.count_nonzero(det & (el == e)))
        row["detached"] = int(det.sum())
        row["fragments"] = int(np.unique(lab[det]).size)
        rows.append(row)
    df = pd.DataFrame(rows)
    if (df["followed"] != len(follow)).any():
        print(f"WARNING: {len(follow)} atoms followed, fewer found in some frames (lost atoms?)")
    df.to_csv(args.output + ".csv", index=False)

    in_frag = np.isin(lab, np.unique(lab[det])) & (lab >= 0)   # with --connect, connecting atoms too
    frag = Counter(O.formulas(np.where(in_frag, lab, -1), el).values())
    pd.DataFrame(sorted(frag.items(), key=lambda kv: -kv[1]),
                 columns=["formula", "count"]).to_csv(args.output + "_fragments.csv", index=False)
    last = rows[-1]
    print("last frame: " + "  ".join(f"{e} {last[f'detached {e}']}" for e in tracked_els)
          + f"  (total {last['detached']}, {last['fragments']} fragments)")
    if frag:
        print("fragments: " + ", ".join(f"{n} x {fm}" for fm, n in frag.most_common(12))
              + (" ..." if len(frag) > 12 else ""))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(6.5, 5.5), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    for k, e in enumerate(tracked_els[:len(COLORS)]):
        O.plot_series(ax, df["time"], df[f"detached {e}"], COLORS[k], e, args.rolling)
    if len(tracked_els) > 1:
        O.plot_series(ax, df["time"], df["detached"], "k", "total", args.rolling, ls="--")
    O.plot_series(ax2, df["time"], df["fragments"], "0.3", None, args.rolling)
    for a in (ax, ax2):
        O.stage_lines(a, df["time"], df["stage"])
        a.set_ylim(bottom=0)
    ax.set_ylabel("detached atoms")
    ax.legend(fontsize=8)
    ax2.set(xlabel=O.time_label(args), ylabel="fragments")
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    print(f"-> {args.output}.csv, {args.output}_fragments.csv, {args.output}.png")


if __name__ == "__main__":
    main()
