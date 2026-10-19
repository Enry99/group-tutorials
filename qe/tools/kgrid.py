#!/usr/bin/env python3
"""Choose a Monkhorst-Pack grid from a k-point spacing, or check a grid.

The same spacing gives the same sampling quality in any cell, so use it to
move a converged grid from bulk to slab, from primitive to conventional cell,
or to a supercell.

Usage:
    python kgrid.py k_14.pwo --grid 14 14 14       # spacing of a known grid
    python kgrid.py slab.pwi --kspacing 0.15       # grid for a target spacing
    python kgrid.py POSCAR --kppra 1000            # Materials Project style density
    python kgrid.py slab.pwi --kspacing 0.15 --slab   # n3 = 1 (vacuum along c)

Spacing convention: dk = |b_i| / n_i in 1/Angstrom, with |b_i| including 2 pi
(same as ASE kspacing and VASP KSPACING). Smaller dk = denser grid.
kppra = k-points per reciprocal atom = n1*n2*n3*natoms (Materials Project uses
~1000 for relaxations).
Any file ASE can read works: .pwo, .pwi with ibrav = 0 (ASE cannot read
ibrav != 0 inputs), POSCAR, .cif, .xyz with a cell.
Needs ase and numpy.
"""
import argparse

import numpy as np
from ase.io import read


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--grid", nargs=3, type=int)
    group.add_argument("--kspacing", type=float, help="1/Angstrom, 2pi included")
    group.add_argument("--kppra", type=float)
    parser.add_argument("--slab", action="store_true", help="force n3 = 1")
    args = parser.parse_args()

    atoms = read(args.structure)
    bnorm = 2.0 * np.pi * np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    print(f"cell lengths (A)    : {np.round(atoms.cell.lengths(), 4)}")
    print(f"|b_i| (1/A, with 2pi): {np.round(bnorm, 4)}")

    if args.grid:
        grid = np.array(args.grid)
    elif args.kspacing:
        grid = np.maximum(1, np.ceil(bnorm / args.kspacing).astype(int))
    else:
        # same spacing along all directions, total ~ kppra / natoms
        target = args.kppra / len(atoms)
        dk = (np.prod(bnorm) / target) ** (1.0 / 3.0)
        grid = np.maximum(1, np.round(bnorm / dk).astype(int))
    if args.slab:
        grid[2] = 1

    spacing = bnorm / grid
    print(f"grid                : {grid[0]} {grid[1]} {grid[2]}")
    print(f"spacing (1/A)       : {np.round(spacing, 4)}"
          + ("  (n3 = 1: ignore the 3rd value)" if args.slab else ""))
    print(f"kppra               : {int(np.prod(grid)) * len(atoms)}")


if __name__ == "__main__":
    main()
