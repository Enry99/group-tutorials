#!/usr/bin/env python3
"""Stack bottom slab / middle (molecules, film; optional) / top slab or tip
into a LAMMPS data file, and write the LAMMPS group definitions and the
force per atom for a normal load.

Usage:
    python interface.py Fe110_ortho.extxyz --flip-top --middle lubricant.extxyz \\
        --load 1.0 --specorder Fe C H O -o Fe_lub_Fe.data
    python interface.py Fe110.extxyz --no-top --middle water.extxyz --specorder Fe O H \\
        -o Fe_water.data            # molecules on a surface (no top slab)
    python interface.py SiN_slab.extxyz --top C_tip.extxyz --isolated-top \\
        --fixed 5 --thermo 10 --top-fixed 2 --top-slide 2 --top-thermo 6 \\
        --load 5.0 --load-area 2124 --specorder Si C O H N -o tip_on_SiN.data

Layout (z grows upwards), group names used in the lammps/ templates:
    layer_top    outermost --top-fixed A of the top: load is applied here
    slide_top    next --top-slide A (optional): moved for sliding
    thermo_top   next --top-thermo A: thermostatted
    surf_top     rest of the top slab, near the interface: NVE
    molecules    the middle part, if any: NVE
    surf_bot     bottom slab near the interface: NVE
    thermo_bot   next --thermo A: thermostatted
    layer_bot    lowest --fixed A of the bottom slab: kept fixed
The thicknesses are measured from the outer surface of each slab. Atoms are
written bottom, middle, top, each sorted along z, so every group is one
contiguous id range ("id <> first last").

What it does (based on create_interf*.py by E. Pedretti):
  1. reads the slabs; the top is a copy of the bottom if not given
     (--no-top: none, e.g. molecules on a surface);
     --flip-top turns it upside down (180 deg about x), so that the same
     surface faces the interface (rectangular or hexagonal cells only);
  2. checks that the lateral cells match (--max-strain, default 0.5 %) and
     strains the top onto the bottom; --isolated-top skips the check and
     centres the top (a tip, a cluster) above the bottom slab, checking that
     it does not touch its periodic images;
  3. puts the middle --gap A above the bottom and the top --gap A above it,
     adds --vacuum A above the top (cell c);
  4. writes the data file (atom types in --specorder: the order of
     pair_coeff), <output>.groups.lmp with the group commands, and prints
     the force per atom of layer_top for --load GPa over the cell area (or
     --load-area A^2, e.g. the contact area of a tip).
Include the groups in the LAMMPS input with: include <output>.groups.lmp
Needs ase, numpy.
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import force_per_atom, sort_z, write_structure  # noqa: E402


def zones(z, edge, widths, from_top):
    """Assign each atom to the first zone (outermost first) whose cumulative
    thickness from the outer edge contains it; the last zone takes the rest."""
    depth = (edge - z) if from_top else (z - edge)
    label = np.full(len(z), len(widths), dtype=int)
    limit = 0.0
    for k, w in enumerate(widths):
        limit += w
        label[(depth < limit) & (label == len(widths))] = k
    return label


def mirror_ok(cell):
    """True if the in-plane lattice maps onto itself under y -> -y, as the
    180 deg rotation about x of --flip-top requires (rectangular cells, and
    hexagonal ones with a1 along x)."""
    m = cell[:2, :2].T                                  # columns: a1, a2 (x, y)
    for v in cell[:2, :2]:
        n = np.linalg.solve(m, [v[0], -v[1]])
        if not np.allclose(n, np.round(n), atol=1e-6):
            return False
    return True


def id_range(mask, offset):
    idx = np.where(mask)[0] + offset + 1
    if len(idx) == 0:
        return None
    if np.any(np.diff(idx) != 1):
        sys.exit("internal error: group is not contiguous")
    return idx[0], idx[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bottom")
    parser.add_argument("--top", help="top slab or tip (default: copy of the bottom)")
    parser.add_argument("--flip-top", action="store_true")
    parser.add_argument("--no-top", action="store_true",
                        help="bottom slab and middle only (e.g. molecules on a surface)")
    parser.add_argument("--isolated-top", action="store_true",
                        help="top is not periodic in xy (tip): no cell check, centred")
    parser.add_argument("--middle", help="molecules or film between the slabs")
    parser.add_argument("--gap", type=float, default=2.5, help="A between the parts")
    parser.add_argument("--vacuum", type=float, default=10.0, help="A above the top")
    parser.add_argument("--fixed", type=float, default=1.5, help="layer_bot thickness (A)")
    parser.add_argument("--thermo", type=float, default=6.5, help="thermo_bot thickness (A)")
    parser.add_argument("--top-fixed", type=float, help="layer_top (default: --fixed)")
    parser.add_argument("--top-slide", type=float, default=0.0, help="slide_top (A), 0 = none")
    parser.add_argument("--top-thermo", type=float, help="thermo_top (default: --thermo)")
    parser.add_argument("--max-strain", type=float, default=0.005)
    parser.add_argument("--load", type=float, help="normal load (GPa)")
    parser.add_argument("--load-area", type=float, help="area (A^2), default: cell area")
    parser.add_argument("--specorder", nargs="+", required=True,
                        help="element order of the atom types (as in pair_coeff)")
    parser.add_argument("-o", "--output", default="interface.data")
    args = parser.parse_args()
    top_fixed = args.fixed if args.top_fixed is None else args.top_fixed
    top_thermo = args.thermo if args.top_thermo is None else args.top_thermo

    bottom = sort_z(read(args.bottom))
    top = None if args.no_top else (read(args.top) if args.top else bottom.copy())
    if args.flip_top and top is not None:
        if not args.isolated_top and not mirror_ok(top.cell[:]):
            sys.exit("--flip-top: the in-plane cell is not symmetric under y -> -y "
                     "(oblique cell); make it orthorhombic first (ortho_cell.py)")
        top.rotate(180, "x", center="COM")
    middle = sort_z(read(args.middle)) if args.middle else None

    cell = bottom.cell[:].copy()
    if np.abs(bottom.cell.angles()[:2] - 90).max() > 1e-4:
        sys.exit("bottom: c must be perpendicular to the surface (slab.py does it)")
    if top is None:
        pass
    elif args.isolated_top:
        shift = bottom.positions[:, :2].mean(0) - top.positions[:, :2].mean(0)
        top.positions[:, :2] += shift
        clear = np.diag(cell)[:2] - np.ptp(top.positions[:, :2], axis=0)
        if clear.min() < 2.0:
            sys.exit(f"the top ({np.ptp(top.positions[:, :2], axis=0).round(1)} A wide) "
                     f"does not fit in the lateral cell {np.diag(cell)[:2].round(1)} A")
        if clear.min() < 10.0:
            print(f"WARNING: the top is only {clear.min():.1f} A from its periodic image")
    else:
        strain = np.abs(top.cell[:2, :2] - cell[:2, :2]).max() / cell[:2, :2].max()
        if strain > args.max_strain:
            sys.exit(f"lateral cells differ by {strain:.2%} (> {args.max_strain:.2%}):\n"
                     f"bottom {np.round(cell[:2, :2], 4).tolist()}\n"
                     f"top    {np.round(top.cell[:2, :2], 4).tolist()}")
        tcell = top.cell[:].copy()
        tcell[:2, :2] = cell[:2, :2]
        top.set_cell(tcell, scale_atoms=True)
        print(f"top strained onto the bottom cell ({strain:.3%})")
    if top is not None:
        top = sort_z(top)

    # stack along z
    bottom.positions[:, 2] -= bottom.positions[:, 2].min()
    zmax = bottom.positions[:, 2].max()
    if middle is not None:
        middle.positions[:, 2] += zmax + args.gap - middle.positions[:, 2].min()
        zmax = middle.positions[:, 2].max()
    # zones from the outer surfaces
    zb = zones(bottom.positions[:, 2], bottom.positions[:, 2].min(),
               [args.fixed, args.thermo], from_top=False)        # 0 layer, 1 thermo, 2 surf
    names_top = []
    if top is not None:
        top.positions[:, 2] += zmax + args.gap - top.positions[:, 2].min()
        widths_top = [top_fixed] + ([args.top_slide] if args.top_slide > 0 else []) + \
            [top_thermo]
        zt = zones(top.positions[:, 2], top.positions[:, 2].max(), widths_top, from_top=True)
        names_top = ["layer_top"] + (["slide_top"] if args.top_slide > 0 else []) + \
            ["thermo_top", "surf_top"]

    parts = [bottom] + ([middle] if middle is not None else []) + ([top] if top else [])
    atoms = parts[0].copy()
    for p in parts[1:]:
        atoms += p
    atoms.set_cell(cell)
    atoms.cell[2] = [0.0, 0.0, atoms.positions[:, 2].max() + args.vacuum]
    atoms.pbc = True
    atoms.wrap()

    nb, nm = len(bottom), (len(middle) if middle is not None else 0)
    groups = []
    for k, name in enumerate(["layer_bot", "thermo_bot", "surf_bot"]):
        groups.append((name, id_range(zb == k, 0)))
    if middle is not None:
        groups.append(("molecules", (nb + 1, nb + nm)))
    for k, name in reversed(list(enumerate(names_top))):
        groups.append((name, id_range(zt == k, nb + nm)))
    missing = [n for n, r in groups if r is None]
    if missing:
        sys.exit(f"empty group(s) {missing}: check the thicknesses")

    write_structure(args.output, atoms, "lammps-data", args.specorder)
    stem = os.path.splitext(args.output)[0]
    lines = [f"# groups for {os.path.basename(args.output)}, written by interface.py"]
    lines += [f"group {n:<11s} id <> {a} {b}" for n, (a, b) in groups]
    present = [n for n, _ in groups]
    lines += ["group sub_bot     union layer_bot thermo_bot surf_bot"]
    if top is not None:
        lines += ["group sub_top     union " + " ".join(n for n in names_top if n in present)]
    inter = ["surf_bot"] + (["molecules"] if middle is not None else []) + \
        (["surf_top"] if top is not None else [])
    lines += ["group intercalate union " + " ".join(inter),
              "group mobile      union thermo_bot intercalate"
              + (" thermo_top" if top is not None else "")
              + (" slide_top" if top is not None and args.top_slide > 0 else "")]

    area = args.load_area or float(np.linalg.norm(np.cross(cell[0], cell[1])))
    if args.load is not None and top is None:
        sys.exit("--load needs a top slab")
    if args.load is not None:
        n_load = groups[[n for n, _ in groups].index("layer_top")][1]
        n_load = n_load[1] - n_load[0] + 1
        f = force_per_atom(args.load, area, n_load)
        lines += [f"variable    load_force equal {-f:.10g}   "
                  f"# eV/A per atom: {args.load} GPa on {area:.2f} A^2, {n_load} atoms"]
    with open(stem + ".groups.lmp", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"{len(atoms)} atoms ({nb} bottom, {nm} middle, {len(top) if top else 0} top), "
          f"cell {np.round(atoms.cell.lengths(), 4)} A -> {args.output}")
    print("\n".join(lines[1:]))
    print(f"-> {stem}.groups.lmp")


if __name__ == "__main__":
    main()
