#!/usr/bin/env python3
"""Turn a cell with an oblique in-plane lattice (hexagonal, bcc(110), ...)
into the smallest orthorhombic supercell, e.g. for sliding MD along x or y.

Usage:
    python ortho_cell.py Fe2O3_hex.extxyz                  # -> Fe2O3_hex_ortho.extxyz
    python ortho_cell.py slab.pwi --repeat 3 2 -f lammps-data --specorder Fe O -o slab.data
    python ortho_cell.py MX-Fe2O3.xyz --square 20          # lateral sides close to 20 A

What it does (based on a script by E. Pedretti):
  1. rotates the structure about z so that a1 lies along x;
  2. finds the integers (n, m), smallest m > 0, with n a1 + m a2 perpendicular
     to a1, and builds the supercell  [a1, n a1 + m a2, a3]
     (hexagonal 120 deg: [1,0,0],[1,2,0],[0,0,1]; 60 deg: [-1,2,0]; bcc(110): [-1,2,0]);
  3. optionally repeats it in-plane (--repeat nx ny, or --square L);
  4. checks that the result is orthogonal and that no atoms overlap;
  5. sorts the atoms along z (group convention; --no-sort to skip).
The number of atoms grows by m (2 for hexagonal). Per-atom data (qe_type,
if_pos, magnetic moments, tags) are carried along.

Requirements: a3 perpendicular to the plane of a1, a2, and a1, a2 in the xy
plane (true for slabs and hexagonal cells from QE, ASE, Materials Project in
the standard setting). Reads .pwi (any ibrav, via qe2xyz.py), .pwo and every
format ASE reads. Writes any ASE format (-f): extxyz (default), lammps-data,
vasp, cif. For a QE input, pass the result to struct_to_qe.py.
Needs ase, numpy.
"""
import argparse
import os
import sys
from fractions import Fraction

import numpy as np
from ase.build import make_supercell
from ase.io import write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe2xyz import min_distance  # noqa: E402
from qe_common import load_structure, sort_z  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--repeat", nargs=2, type=int, default=(1, 1), metavar=("NX", "NY"))
    group.add_argument("--square", type=float, metavar="L",
                       help="repeat so that both lateral sides are close to L (A)")
    parser.add_argument("--mmax", type=int, default=12, help="largest multiplier tried")
    parser.add_argument("-f", "--format", default="extxyz")
    parser.add_argument("--no-sort", action="store_true",
                        help="keep the supercell order instead of sorting along z")
    parser.add_argument("--specorder", nargs="+",
                        help="lammps-data: element order of the atom types (as in pair_coeff)")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    atoms = load_structure(args.structure)
    cell = atoms.cell[:]
    a1, a2, a3 = cell
    if abs(a1[2]) > 1e-6 or abs(a2[2]) > 1e-6:
        sys.exit("a1 and a2 must lie in the xy plane")
    if abs(np.dot(a3, a1)) > 1e-6 or abs(np.dot(a3, a2)) > 1e-6:
        sys.exit("a3 must be perpendicular to a1 and a2 (slab or hexagonal-type cell)")

    # 1. a1 along x
    angle = np.degrees(np.arctan2(a1[1], a1[0]))
    atoms.rotate(-angle, "z", rotate_cell=True)
    a1, a2 = atoms.cell[0], atoms.cell[1]

    # 2. n a1 + m a2 perpendicular to a1:  n = -m (a1.a2)/|a1|^2 must be an integer
    ratio = np.dot(a1, a2) / np.dot(a1, a1)
    for m in range(1, args.mmax + 1):
        n = -m * ratio
        if abs(n - round(n)) < 1e-4:
            n = int(round(n))
            break
    else:
        frac = Fraction(ratio).limit_denominator(100)
        sys.exit(f"no orthogonal supercell with m <= {args.mmax} "
                 f"(a1.a2/|a1|^2 = {ratio:.6f} ~ {frac}); raise --mmax?")
    if a2[1] * m < 0:  # keep the new b along +y (right-handed cell)
        n, m = -n, -m
    P = np.array([[1, 0, 0], [n, m, 0], [0, 0, 1]])
    ortho = make_supercell(atoms, P)

    # 3. in-plane repetition
    lx, ly = ortho.cell.lengths()[:2]
    if args.square:
        nx, ny = max(1, round(args.square / lx)), max(1, round(args.square / ly))
    else:
        nx, ny = args.repeat
    if (nx, ny) != (1, 1):
        ortho = ortho.repeat((nx, ny, 1))
    ortho.wrap()
    if not args.no_sort:
        ortho = sort_z(ortho)          # group convention

    # 4. checks
    angles = ortho.cell.angles()
    if np.abs(angles - 90.0).max() > 1e-4:
        sys.exit(f"result is not orthogonal: angles {angles}")
    if len(ortho) != len(atoms) * abs(m) * nx * ny:
        sys.exit("unexpected number of atoms: duplicated or lost atoms")
    d0, d1 = min_distance(atoms), min_distance(ortho)
    if abs(d0 - d1) > 1e-6:
        sys.exit(f"shortest distance changed ({d0:.4f} -> {d1:.4f} A): something is wrong")

    out = args.output or os.path.splitext(args.structure)[0] + "_ortho." + (
        "extxyz" if args.format == "extxyz" else args.format.replace("-", "_"))
    if args.format == "lammps-data":
        kwargs = {"masses": True, "atom_style": "atomic"}
        if args.specorder:
            kwargs["specorder"] = args.specorder
        write(out, ortho, format=args.format, **kwargs)
    else:
        write(out, ortho, format=args.format)
    lx, ly, lz = ortho.cell.lengths()
    print(f"supercell matrix [[1,0,0],[{n},{m},0],[0,0,1]] x ({nx},{ny},1): "
          f"{len(atoms)} -> {len(ortho)} atoms")
    print(f"orthorhombic cell {lx:.4f} x {ly:.4f} x {lz:.4f} A  -> {out}")
    if args.square is None and (nx, ny) == (1, 1):
        r = ly / lx
        print(f"b/a = {r:.4f}: e.g. --repeat {max(1, round(r * 3))} 3 gives "
              f"{lx * max(1, round(r * 3)):.2f} x {ly * 3:.2f} A")


if __name__ == "__main__":
    main()
