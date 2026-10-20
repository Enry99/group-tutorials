#!/usr/bin/env python3
"""Repeat a bulk cell and scale it to a target mass density (g/cm^3), e.g.
for the starting structure of a melt-quench.

Usage:
    python density.py C_doped.extxyz --target 3.0                 # all three sides
    python density.py TiO2_rutile.cif --repeat 4 3 6 --target 3.8
    python density.py TiO2_rutile.cif --repeat 4 3 6 --xy 19.856 21.061 --target 3.8 \\
        -f lammps-data --specorder Fe Ti O -o TiO2_for_Fe110.data
    python density.py Si3N4.cif --repeat 3 3 6 --target 2.9 --remove Si3N4   # cell kept

What it does:
  1. reads the structure and repeats it (--repeat nx ny nz);
  2. --xy LX LY: sets the lateral sides to those of the slab it will sit on
     (orthorhombic cells only; atoms scaled with the cell); then only c is
     scaled to reach the target density, so the lateral cell is kept;
     without --xy all three sides are scaled by the same factor;
  3. --remove FORMULA: the cell is not scaled (after --xy, if given); the
     density is lowered by removing whole formula units (e.g. Si3N4: 3 Si
     and 4 N each time) at random positions (--seed), so the composition
     does not change (based on the SiN and SiO2 notebooks of the group);
  4. prints the density and the cell before and after, writes the result.
The target density is that of the amorphous phase you want (experiment or
literature), not of the crystal. Strained or compressed starting cells, or
vacancies, are fine for a melt-quench: the melt forgets them. Needs ase,
numpy.
"""
import argparse
import os
import sys

import numpy as np
from ase import Atoms
from ase.io import read

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import density, extension, write_structure  # noqa: E402


def report(tag, atoms):
    a, b, c = atoms.cell.lengths()
    print(f"{tag}: {len(atoms)} atoms, cell {a:.4f} x {b:.4f} x {c:.4f} A, "
          f"density {density(atoms):.3f} g/cm^3")


def remove_units(atoms, formula, target, rng):
    """Delete whole formula units at random positions, as many as bring the
    density closest to target (it can only decrease)."""
    unit = Atoms(formula)
    n = round((1 - target / density(atoms)) * atoms.get_masses().sum()
              / unit.get_masses().sum())
    if n <= 0:
        sys.exit(f"--remove lowers the density: the start ({density(atoms):.3f} g/cm^3) "
                 "is already at or below the target")
    symbols = np.array(atoms.get_chemical_symbols())
    gone = []
    for el, k in unit.symbols.formula.count().items():
        idx = np.flatnonzero(symbols == el)
        if n * k > len(idx):
            sys.exit(f"--remove {formula}: {n} units need {n * k} {el}, there are {len(idx)}")
        gone.extend(rng.choice(idx, n * k, replace=False))
    del atoms[gone]
    print(f"removed {n} units of {formula} ({len(gone)} atoms)")
    return atoms


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure")
    parser.add_argument("--target", type=float, required=True, help="g/cm^3")
    parser.add_argument("--repeat", nargs=3, type=int, default=(1, 1, 1))
    parser.add_argument("--xy", nargs=2, type=float, metavar=("LX", "LY"),
                        help="lateral sides (A); then only c is scaled")
    parser.add_argument("--remove", metavar="FORMULA",
                        help="keep the cell, remove random formula units (e.g. Si3N4)")
    parser.add_argument("--seed", type=int, default=0, help="for --remove")
    parser.add_argument("-f", "--format", default="extxyz")
    parser.add_argument("--specorder", nargs="+",
                        help="lammps-data: element order of the atom types (as in pair_coeff)")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    atoms = read(args.structure) * tuple(args.repeat)
    atoms.pbc = True
    report("start", atoms)
    if args.xy:
        if np.abs(atoms.cell.angles() - 90).max() > 1e-6:
            sys.exit("--xy needs an orthorhombic cell (qe/tools/ortho_cell.py first)")
        atoms.set_cell([args.xy[0], args.xy[1], atoms.cell[2, 2]], scale_atoms=True)
        report("lateral cell set", atoms)
    if args.remove:
        atoms = remove_units(atoms, args.remove, args.target, np.random.default_rng(args.seed))
    else:
        if args.xy:
            cell = atoms.cell[:]
            cell[2] *= density(atoms) / args.target
        else:
            cell = atoms.cell[:] * (density(atoms) / args.target) ** (1 / 3)
        atoms.set_cell(cell, scale_atoms=True)
    atoms.wrap()
    report("final", atoms)

    stem = os.path.splitext(os.path.basename(args.structure))[0]
    ext = extension(args.format)
    out = args.output or f"{stem}_{args.target:.2f}gcm3.{ext}"
    write_structure(out, atoms, args.format, args.specorder)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
