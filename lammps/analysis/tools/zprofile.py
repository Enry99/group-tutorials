#!/usr/bin/env python3
"""Density profiles along z from LAMMPS text dumps, per element or per
group, averaged over a window of frames (no OVITO needed).

Usage:
    python zprofile.py 02_slide.lammpstrj --last 20
    python zprofile.py 01_load.lammpstrj 02_slide.lammpstrj --tmin 60 --dt 0.0005 \\
        --elements Fe O H --mass --sigma 0.5 --shade Fe -o water_film
    python zprofile.py 02_slide.lammpstrj --groups Fe_water_Fe.groups.lmp \\
        --by-group surf_bot molecules surf_top --last 20 --avechunk

What it does (based on density_profile_from_dump.py, E. Pedretti):
  1. reads the frames of the dumps (files in the order given, a frame
     repeated at a restart counted once); the window is the last --last
     frames, or --tmin/--tmax in time (step x --dt per stage, stages laid
     end to end after a reset_timestep), every --every frame;
  2. bins z (--bin, default 0.5 A, from the lowest box bound); z is the
     unwrapped column (zu, else z) folded into the box if z is periodic;
  3. per frame, number density of each element (column 'element', or
     'type' with --type-names) or of each group of a groups file (needs
     'id' in the dump): N(bin) / (lx ly dz), in A^-3, or g/cm^3 with
     --mass; frames averaged with equal weight, as fix ave/chunk with
     'norm sample';
  4. writes <out>.csv (z and one column per element or group) and
     <out>.png (z vertical; --sigma smooths the plot only, in A; --shade X
     greys the regions where element X, averaged over 1.5 A, is above a
     quarter of its maximum, e.g. the slabs); --avechunk also writes dens_profile_<X> in the format of
     fix ave/chunk (density/number).
Needs numpy, pandas, matplotlib, ase (masses, colours), scipy (--sigma).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp_io as L  # noqa: E402

AMU_A3_TO_G_CM3 = 1.66053906660    # 1 amu/A^3 in g/cm^3
COLORS = {"H": (0.85, 0.40, 0.70), "C": (0.20, 0.20, 0.25), "O": (0.70, 0.0, 0.0),
          "S": "gold", "Fe": (0.40, 0.40, 1.0), "Ti": "teal"}     # his choices; jmol otherwise


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dumps", nargs="+")
    parser.add_argument("--last", type=int, help="last N frames")
    parser.add_argument("--tmin", type=float)
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--dt", nargs="+", type=float, help="timestep per stage (for --tmin/--tmax)")
    parser.add_argument("--every", type=int, default=1)
    parser.add_argument("--bin", type=float, default=0.5, help="A")
    parser.add_argument("--elements", nargs="+", help="default: all in the dump")
    parser.add_argument("--type-names", nargs="+", help="element of type 1, 2, ... (no 'element' column)")
    parser.add_argument("--groups", help="groups file (structures/tools/interface.py)")
    parser.add_argument("--by-group", nargs="+", metavar="GROUP")
    parser.add_argument("--mass", action="store_true", help="g/cm^3 instead of A^-3")
    parser.add_argument("--sigma", type=float, default=0.0, help="Gaussian smoothing of the plot (A)")
    parser.add_argument("--shade", help="element whose dense regions are greyed (slabs)")
    parser.add_argument("--avechunk", action="store_true", help="also dens_profile_<X> files")
    parser.add_argument("-o", "--output", default="zprofile")
    args = parser.parse_args()

    idx = L.index_dump(L.expand(args.dumps))
    steps = np.array([s for _, _, s in idx])
    sel, _ = L.window(steps, args.dt, args.tmin, args.tmax, args.last, args.every)
    print(f"{len(idx)} frames in the dumps; window: {len(sel)} frames, steps "
          f"{steps[sel[0]]}-{steps[sel[-1]]}")

    from ase.data import atomic_masses, atomic_numbers
    group_ids = None
    if args.by_group:
        if not args.groups:
            sys.exit("--by-group needs --groups FILE")
        group_ids = {g: L.ids_of(g, args.groups) for g in args.by_group}

    counts, ncounts, areas, zgrid = {}, {}, [], None
    for k in sel:
        fn, pos, _ = idx[k]
        f = L.read_frame(fn, pos)
        d = f["cols"]
        zc = "zu" if "zu" in d else "z"
        z = d[zc].to_numpy(dtype=float)
        lz = f["hi"][2] - f["lo"][2]
        if f["pbc"][2]:
            z = f["lo"][2] + np.mod(z - f["lo"][2], lz)
        if zgrid is None:                                # fixed grid from the first frame of the window
            z0 = f["lo"][2]
            nb = int(np.ceil(lz / args.bin)) + 1
            zgrid = z0 + args.bin * (np.arange(nb) + 0.5)
        ib = np.clip(np.floor((z - z0) / args.bin).astype(int), 0, len(zgrid) - 1)
        areas.append((f["hi"][0] - f["lo"][0]) * (f["hi"][1] - f["lo"][1]))
        if group_ids is not None:
            if "id" not in d:
                sys.exit("--by-group needs 'id' in the dump")
            ids = d["id"].to_numpy()
            masks = {g: np.isin(ids, v) for g, v in group_ids.items()}
            el = (d["element"].to_numpy() if "element" in d else None)
        else:
            if "element" in d:
                el = d["element"].astype(str).to_numpy()
            elif args.type_names:
                el = np.array(args.type_names)[d["type"].to_numpy() - 1]
            else:
                sys.exit("no 'element' column: give --type-names")
            wanted = args.elements or sorted(set(el), key=lambda e: atomic_numbers.get(e, 0))
            masks = {e: el == e for e in wanted}
        for name, m in masks.items():
            if args.mass:
                if el is None:
                    sys.exit("--mass by group needs the 'element' column")
                w = np.array([atomic_masses[atomic_numbers[e]] for e in el[m]])
            else:
                w = None
            c = np.bincount(ib[m], weights=w, minlength=len(zgrid))
            counts.setdefault(name, []).append(c / (areas[-1] * args.bin))
            if args.avechunk:                               # number density, also with --mass
                nc = np.bincount(ib[m], minlength=len(zgrid))
                ncounts.setdefault(name, []).append(nc / (areas[-1] * args.bin))

    unit = "g/cm^3" if args.mass else "A^-3"
    conv = AMU_A3_TO_G_CM3 if args.mass else 1.0
    prof = {name: np.sum(c, axis=0) / len(sel) * conv   # frames without the element count as 0
            for name, c in counts.items()}
    out = pd.DataFrame({"z": zgrid, **{f"{n} ({unit})": p for n, p in prof.items()}})
    out.to_csv(args.output + ".csv", index=False)
    for n, p in prof.items():
        print(f"  {n:>12}: max {p.max():.4g} {unit}" + ("   <- empty" if p.max() == 0 else ""))
    if args.avechunk:
        for n, c in ncounts.items():
            dens = np.sum(c, axis=0) / len(sel)
            mean_n = dens * np.mean(areas) * args.bin
            with open(f"dens_profile_{n}", "w") as fh:
                fh.write(f"# Chunk-averaged data for {n} (zprofile.py)\n")
                fh.write("# Timestep Number-of-chunks Total-count\n")
                fh.write("# Chunk Coord1 Ncount density/number\n")
                fh.write(f"{steps[sel[-1]]} {len(zgrid)} {mean_n.sum():g}\n")
                for j in range(len(zgrid)):
                    fh.write(f"  {j + 1} {zgrid[j]:g} {mean_n[j]:g} {dens[j]:g}\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase.data.colors import jmol_colors
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(figsize=(4.5, 7))
    styles = ["-", "--", "-.", ":"]
    for k, (n, p) in enumerate(prof.items()):
        if args.sigma:
            from scipy.ndimage import gaussian_filter1d
            p = gaussian_filter1d(p, args.sigma / args.bin, mode="nearest")
        col = COLORS.get(n, jmol_colors[atomic_numbers[n]] if n in atomic_numbers else None)
        ax.plot(p, zgrid, ls=styles[k % 4], lw=2, color=col, label=n)
    if args.shade and args.shade in prof:              # atomic layers averaged out first
        from scipy.ndimage import gaussian_filter1d
        sm = gaussian_filter1d(prof[args.shade], 1.5 / args.bin, mode="nearest")
        dense = sm > 0.25 * sm.max()
        edges = np.flatnonzero(np.diff(np.r_[0, dense.astype(int), 0]))
        for a, b in zip(edges[::2], edges[1::2]):
            ax.axhspan(zgrid[a] - args.bin / 2, zgrid[b - 1] + args.bin / 2,
                       color="0.85", zorder=0, lw=0)
        others = [p[~dense].max() for n, p in prof.items() if n != args.shade]
        if others and max(others) > 0:                 # x range set by the film, slabs clipped
            ax.set_xlim(0, 1.1 * max(others))
    ax.set(xlabel=f"density ({unit})", ylabel="z (Å)", ylim=(zgrid[0], zgrid[-1]))
    ax.set_xlim(left=0)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    print(f"-> {args.output}.csv, {args.output}.png")


if __name__ == "__main__":
    main()
