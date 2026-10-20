#!/usr/bin/env python3
"""Convert a normal pressure (GPa) into the force per atom (eV/A) for
LAMMPS fix addforce, and back.

Usage:
    python load_force.py 1.0 --cell interface.data --natoms 72    # area of the cell
    python load_force.py 5.0 --radius 15 --natoms 317             # flat face of a tip
    python load_force.py --force -0.0362516 --cell interface.data --natoms 72

F_atom = P A / N, with A the area on which the pressure acts (the xy area of
the cell for a periodic slab, pi R^2 or --area for a tip) and N the atoms of
the group the force is applied to. 1 eV/A = 1.602176634e-9 N. Needs ase.
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import EV_PER_A_TO_N, force_per_atom  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pressure", type=float, nargs="?", help="GPa")
    parser.add_argument("--force", type=float, help="eV/A per atom: print the pressure")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cell", help="structure whose xy cell area is used")
    group.add_argument("--area", type=float, help="A^2")
    group.add_argument("--radius", type=float, help="A: area pi R^2")
    parser.add_argument("--natoms", type=int, required=True)
    args = parser.parse_args()

    if args.cell:
        kw = {"format": "lammps-data"} if args.cell.endswith(".data") else {}
        cell = read(args.cell, **kw).cell
        area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    else:
        area = args.area if args.area else np.pi * args.radius ** 2
    if args.force is not None:
        p = abs(args.force) * args.natoms * EV_PER_A_TO_N / (area * 1e-20) / 1e9
        print(f"{args.force} eV/A on {args.natoms} atoms over {area:.2f} A^2 = {p:.4f} GPa")
    elif args.pressure is not None:
        f = force_per_atom(args.pressure, area, args.natoms)
        print(f"{args.pressure} GPa over {area:.2f} A^2 on {args.natoms} atoms: "
              f"{f:.10g} eV/A per atom (use -{f:.10g} to push down)")
    else:
        sys.exit("give a pressure or --force")


if __name__ == "__main__":
    main()
