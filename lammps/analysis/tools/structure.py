#!/usr/bin/env python3
"""Local structure from LAMMPS text dumps (OVITO): partial RDFs with a check
of the bond cutoffs, coordination numbers, bond-angle distributions.

Usage:
    python structure.py meltquench.lammpstrj
    python structure.py meltquench.lammpstrj --last 10 --angles C-C-C -o aC
    python structure.py anneal.lammpstrj --zmin 12 --zmax 30 --angles N-Si-N Si-N-Si \\
        --cutoff Si-N:2.37

What it does:
  1. frames: the last one by default, or a window (--last N, --tmin/--tmax
     with --dt, --every); quantities averaged over the frames;
  2. partial RDFs g_AB(r) up to --rmax (OVITO CoordinationAnalysisModifier,
     partial=True: each partial tends to 1); for each pair, the first peak
     and the valley after it (first minimum as in rdf_check.py, E.
     Pedretti, on g smoothed over 0.05 A, widened to where g stays within
     5 % of the peak height of the minimum), the cutoff in use and
     g(cutoff)/g(peak) and the A-B bonds per frame; a cutoff outside the
     valley is flagged: below the first peak with no bonds it means no A-B
     bonds (right for a non-bonded contact, e.g. H-H in water), with bonds
     a small peak below the first one (look at the plot), on a peak
     it cuts the first shell, beyond the valley it takes in the next
     shell (e.g. Fe-Fe from the covalent radii in bcc Fe);
     The figure shows the cutoff in use (green dashed), the valley (red
     band) and, where cutoffs.py has two values, the other one (purple);
  3. bonds with the pairwise cutoffs of cutoffs.py (--cutoff A-B:R
     overrides; --alternative for the ovito_cho.py set); coordination of
     each centre atom: distribution of the total number of bonded
     neighbours per element (a-C: CN 2, 3, 4 = sp, sp2, sp3), mean
     number of B neighbours of A for every pair; centres restricted with
     --zmin/--zmax (A, wrapped z) or --centres (ids, groups with --groups,
     elements);
  4. --angles A-B-C: angle at B between bonds B-A and B-C, histogram at
     0.1 deg, plotted smoothed by a Gaussian of --angle-sigma deg;
  5. writes <out>_rdf.csv/png, <out>_cn.csv/png and <out>_angles.csv/png.
Needs ovito, numpy, pandas, matplotlib, scipy.
"""
import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cutoffs as C  # noqa: E402
import lmp_io as L  # noqa: E402
import ovito_io as O  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def valley(r, g):
    """First peak of g(r) and the valley after it: g smoothed over 0.05 A,
    first maximum above 30 % of the highest one, first minimum after it, and the range around it where g is
    within 5 % of the peak height of that minimum. Returns (peak, lo, hi);
    lo, hi are None if g does not come down after the peak."""
    from scipy.ndimage import gaussian_filter1d
    if g.max() <= 0:
        return None, None, None
    gs = gaussian_filter1d(g, 0.05 / (r[1] - r[0]), mode="nearest")
    top = np.flatnonzero((gs[1:-1] >= gs[:-2]) & (gs[1:-1] > gs[2:])
                         & (gs[1:-1] > 0.3 * gs.max())) + 1
    k = int(top[0]) if len(top) else int(np.argmax(gs))     # first peak above 30 % of the highest
    for t in range(k + 1, len(gs) - 1):
        if gs[t] <= gs[t + 1] and gs[t] < 0.4 * gs[k]:
            low = gs <= gs[t] + 0.05 * gs[k]
            a = b = t
            while a > k and low[a - 1]:
                a -= 1
            while b < len(gs) - 1 and low[b + 1]:
                b += 1
            return r[k], r[a], r[b]
    return r[k], None, None


def angles_at(i, j, vec, el, centres, triplets):
    """Bond angles (deg) for each triplet A-B-C, B a centre atom."""
    out = {t: [] for t in triplets}
    ci = np.r_[i, j]                                     # half-bonds: centre, neighbour, vector
    ni = np.r_[j, i]
    vv = np.r_[vec, -vec]
    m = centres[ci]
    ci, ni, vv = ci[m], ni[m], vv[m]
    order = np.argsort(ci, kind="stable")
    ci, ni, vv = ci[order], ni[order], vv[order]
    starts = np.flatnonzero(np.r_[True, ci[1:] != ci[:-1]])
    deg = np.diff(np.r_[starts, len(ci)])
    for d in np.unique(deg):
        if d < 2:
            continue
        rows = starts[deg == d][:, None] + np.arange(d)   # (centres, d)
        for p, q in itertools.combinations(range(d), 2):
            a, b = rows[:, p], rows[:, q]
            cos = np.einsum("ij,ij->i", vv[a], vv[b]) / (
                np.linalg.norm(vv[a], axis=1) * np.linalg.norm(vv[b], axis=1))
            ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
            ea, eb, ec = el[ni[a]], el[ci[a]], el[ni[b]]
            for t in triplets:
                A, B, Cc = t
                hit = (eb == B) & (((ea == A) & (ec == Cc)) | ((ea == Cc) & (ec == A)))
                out[t].append(ang[hit])
    return {t: np.concatenate(v) if v else np.zeros(0) for t, v in out.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    O.add_common_args(parser, last_default=1)
    parser.add_argument("--rmax", type=float, default=5.0, help="RDF range (A)")
    parser.add_argument("--rbins", type=int, default=250)
    parser.add_argument("--zmin", type=float, help="centres with z above (A)")
    parser.add_argument("--zmax", type=float)
    parser.add_argument("--centres", help="centre atoms: ids 'lo:hi', 'a,b', groups 'g1+g2', "
                        "or elements 'C+H'")
    parser.add_argument("--groups", help="groups file for --centres")
    parser.add_argument("--angles", nargs="+", default=[], metavar="A-B-C")
    parser.add_argument("--angle-sigma", type=float, default=1.0, help="deg, plot only")
    parser.add_argument("-o", "--output", default="structure")
    args = parser.parse_args()

    pipeline, steps, els = O.load(args.dumps, args.type_names)
    frames, t, _ = O.window(steps, args)
    print(f"{len(frames)} frame(s): steps {steps[frames[0]]}-{steps[frames[-1]]}")
    triplets = [tuple(a.split("-")) for a in args.angles]
    for tr in triplets:
        if len(tr) != 3 or not set(tr) <= set(els):
            sys.exit(f"--angles {'-'.join(tr)}: three elements of {els}")

    from ovito.modifiers import CoordinationAnalysisModifier
    pipeline.modifiers.append(CoordinationAnalysisModifier(
        cutoff=args.rmax, number_of_bins=args.rbins, partial=True))
    table = O.add_bonds(pipeline, els, args)

    rdf, cn_tot, cn_pair, n_centres, n_bonds = None, {}, {}, {}, {}
    ang = {tr: [] for tr in triplets}
    for k in frames:
        data = pipeline.compute(k)
        xy = data.tables["coordination-rdf"].xy()
        names = list(data.tables["coordination-rdf"].y.component_names)
        rdf = xy if rdf is None else rdf + xy
        el = O.elements_of(data)
        n = len(el)
        centres = np.ones(n, dtype=bool)
        z = np.asarray(data.particles.positions)[:, 2]
        if args.zmin is not None:
            centres &= z >= args.zmin
        if args.zmax is not None:
            centres &= z <= args.zmax
        if args.centres:
            if set(args.centres.split("+")) <= set(els):
                centres &= np.isin(el, args.centres.split("+"))
            else:
                centres &= np.isin(O.ids_of_frame(data), L.ids_of(args.centres, args.groups))
        i, j, vec = O.bond_list(data)
        for (a, b), _ in table.items():
            n_bonds[(a, b)] = n_bonds.get((a, b), 0) + np.count_nonzero(
                ((el[i] == a) & (el[j] == b)) | ((el[i] == b) & (el[j] == a)))
        for a in els:
            ca = centres & (el == a)
            n_centres[a] = n_centres.get(a, 0) + int(ca.sum())
            if not ca.any():
                continue
            tot = np.bincount(i[ca[i]], minlength=n) + np.bincount(j[ca[j]], minlength=n)
            cn_tot[a] = cn_tot.get(a, 0) + np.bincount(np.minimum(tot[ca], 24), minlength=25)
            for b in els:
                nb = (np.count_nonzero(ca[i] & (el[j] == b)) + np.count_nonzero(ca[j] & (el[i] == b)))
                cn_pair[(a, b)] = cn_pair.get((a, b), 0) + nb
        for tr, v in angles_at(i, j, vec, el, centres, triplets).items():
            ang[tr].append(v)
    rdf /= len(frames)

    # ---- RDF and cutoff check
    r = rdf[:, 0]
    print(f"\n{'pair':>7} {'peak (A)':>9} {'valley (A)':>11} {'cutoff':>7} {'g(rc)/g(pk)':>11} "
          f"{'bonds/frame':>11}")
    hints = {}
    for k, name in enumerate(names):
        a, b = name.split("-")
        g = rdf[:, k + 1]
        peak, lo, hi = valley(r, g)
        rc, src = C.cutoff(a, b, C.parse(args.cutoff), args.alternative)
        if peak is None:
            print(f"{name:>7}   no pairs within {args.rmax} A")
            continue
        ratio = np.interp(rc, r, g) / g.max()
        nb = n_bonds.get((a, b), n_bonds.get((b, a), 0)) / len(frames)
        if lo is not None and lo <= rc <= hi:
            flag = ""
        elif rc < peak and ratio < 0.05 and nb == 0:
            flag = "   no bonds: right if A and B do not bond"
        elif rc < peak and ratio < 0.05:
            flag = "   <- check: bonds below the first peak (a small peak before it?)"
        elif rc < (hi if hi is not None else np.inf):
            flag = "   <- check: the cutoff cuts the first peak"
        else:
            flag = "   <- check: beyond the valley, takes in the next shell"
        vtxt = f"{lo:.2f}-{hi:.2f}" if lo is not None else "n/a"
        print(f"{name:>7} {peak:>9.2f} {vtxt:>11} {rc:>7.2f} {ratio:>11.2f} {nb:>11.1f}"
              f"   ({src}){flag}")
        hints[name] = (lo, hi, rc)
    pd.DataFrame(rdf, columns=["r"] + names).to_csv(args.output + "_rdf.csv", index=False)

    # ---- coordination
    rows = []
    print(f"\ncoordination (bonds within the cutoffs; {len(frames)} frame(s))")
    for a in els:
        if a not in cn_tot:
            continue
        h = cn_tot[a] / cn_tot[a].sum()
        mean = (np.arange(len(h)) * h).sum()
        shown = "  ".join(f"{c}: {100 * f:.1f}%" for c, f in enumerate(h) if f >= 0.001)
        print(f"  {a:>2} ({n_centres[a] // len(frames)} centres): mean {mean:.2f}   {shown}")
        per = "  ".join(f"{b} {cn_pair[(a, b)] / n_centres[a]:.2f}" for b in els)
        print(f"       neighbours per {a}: {per}")
        for c, f in enumerate(h):
            rows.append({"centre": a, "neighbour": "all", "n": c, "fraction": f})
        for b in els:
            rows.append({"centre": a, "neighbour": b, "n": "mean",
                         "fraction": cn_pair[(a, b)] / n_centres[a]})
    pd.DataFrame(rows).to_csv(args.output + "_cn.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    nc = min(3, len(names))
    nr = int(np.ceil(len(names) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(3.6 * nc, 2.6 * nr), squeeze=False)
    for k, name in enumerate(names):
        ax = axes.flat[k]
        ax.plot(r, rdf[:, k + 1], color=COLORS[0], lw=1.2)
        ax.axhline(1, color="0.7", lw=0.6)
        if name in hints:
            lo, hi, rc = hints[name]
            ax.axvline(rc, color="tab:green", ls="--", lw=1)
            if lo is not None:
                ax.axvspan(lo, hi, color="tab:red", alpha=0.12, lw=0)
            a, b = name.split("-")
            alt = C.ALTERNATIVE.get((a, b), C.ALTERNATIVE.get((b, a)))
            other = alt if not args.alternative else C.TABLE.get((a, b), C.TABLE.get((b, a)))
            if other is not None and other != rc:
                ax.axvline(other, color="tab:purple", ls="-.", lw=1)
        ax.set_title(name, fontsize=10)
        ax.set_xlim(0, args.rmax)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    fig.supxlabel("r (Å)")
    fig.supylabel("g(r)")
    fig.tight_layout()
    fig.savefig(args.output + "_rdf.png", dpi=200)
    outs = [f"{args.output}_rdf.csv/png", f"{args.output}_cn.csv"]

    present = [a for a in els if a in cn_tot]
    fig, axes = plt.subplots(1, len(present), figsize=(3.2 * len(present), 3), squeeze=False)
    for ax, a in zip(axes[0], present):
        h = cn_tot[a] / cn_tot[a].sum()
        nz = np.flatnonzero(h >= 0.001)
        cs = np.arange(nz.min(), nz.max() + 1)
        ax.bar(cs, 100 * h[cs], color=COLORS[0], width=0.7)
        ax.set(title=f"{a}", xlabel="bonded neighbours", xticks=cs)
    axes[0, 0].set_ylabel("atoms (%)")
    fig.tight_layout()
    fig.savefig(args.output + "_cn.png", dpi=200)
    outs[-1] += "/png"

    if triplets:
        from scipy.ndimage import gaussian_filter1d
        edges = np.arange(0, 180.0001, 0.1)
        mid = 0.5 * (edges[1:] + edges[:-1])
        out = {"angle_deg": mid}
        fig, ax = plt.subplots(figsize=(5.5, 4))
        print("\nbond angles")
        for k, tr in enumerate(triplets):
            v = np.concatenate(ang[tr])
            name = "-".join(tr)
            if len(v) == 0:
                print(f"  {name}: none")
                continue
            h, _ = np.histogram(v, edges, density=True)
            out[name] = h
            sm = gaussian_filter1d(h, args.angle_sigma / 0.1)
            out[name + " smoothed"] = sm
            ax.plot(mid, sm, color=COLORS[k % 4], lw=1.6, label=f"{name} ({len(v) // len(frames)} per frame)")
            print(f"  {name}: {len(v) // len(frames)} per frame, mean {v.mean():.1f} deg, "
                  f"most frequent {mid[np.argmax(sm)]:.1f} deg")
        ax.set(xlabel="angle (deg)", ylabel="probability density (1/deg)", xlim=(40, 180))
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.output + "_angles.png", dpi=200)
        pd.DataFrame(out).to_csv(args.output + "_angles.csv", index=False)
        outs.append(f"{args.output}_angles.csv/png")
    print("-> " + ", ".join(outs))


if __name__ == "__main__":
    main()
