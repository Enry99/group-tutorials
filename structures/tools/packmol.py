#!/usr/bin/env python3
"""Fill a box with molecules at a target density using packmol, e.g. a
lubricant or water film between two slabs, or a liquid around a flake.

Usage:
    python packmol.py H2O --box 19.856 21.061 20 --density 1.0
    python packmol.py C8H18.xyz --cell-from slab.extxyz --thickness 25 --density 0.70
    python packmol.py H2O O2:10 --box 22 22 22 --density 1.0 --solute flake.extxyz

Molecules: a name of the ASE g2 set (H2O, CH4, NH3, O2, H2, CO2, ...) or a
file (xyz, extxyz, ...). The first molecule fills the box at --density
(g/cm^3); the others are added in fixed numbers (NAME:N).

What it does (based on write_packmol.py by E. Pedretti):
  1. box: --box LX LY LZ, or the lateral cell of a slab (--cell-from,
     orthorhombic) and --thickness;
  2. number of molecules = density x free volume / molecular mass; with
     --solute, the solute (kept fixed at the box centre) and its bounding
     box plus 1.5 A are subtracted from the free volume;
  3. writes packmol.inp, runs packmol (on PATH or --packmol) and writes
     <output>.extxyz with the cell of the box. The molecules are packed
     --tolerance/2 away from the faces, so that periodic images do not
     overlap.
Look at the result: packmol prints "Success!" only if the tolerance
(default 2.0 A) is satisfied everywhere. Needs ase, numpy and packmol
(https://github.com/m3g/packmol).
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from ase.build import molecule
from ase.io import read, write

AMU_TO_G = 1.66053906660e-24


def load_molecule(spec):
    name = spec.split(":")[0]
    if os.path.exists(name):
        return os.path.splitext(os.path.basename(name))[0], read(name)
    return name, molecule(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("molecules", nargs="+", help="first: fills the box; others NAME:N")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--box", nargs=3, type=float, metavar=("LX", "LY", "LZ"))
    group.add_argument("--cell-from", help="slab whose lateral (orthorhombic) cell is used")
    parser.add_argument("--thickness", type=float, help="LZ with --cell-from")
    parser.add_argument("--density", type=float, required=True, help="g/cm^3")
    parser.add_argument("--solute", help="structure kept fixed at the box centre")
    parser.add_argument("--tolerance", type=float, default=2.0)
    parser.add_argument("--seed", type=int,
                        help="packmol seed (default: packmol's own, reproducible; -1: random)")
    parser.add_argument("--packmol", default="packmol")
    parser.add_argument("-o", "--output", default="molecules")
    args = parser.parse_args()

    if args.cell_from:
        slab = read(args.cell_from)
        if np.abs(slab.cell.angles() - 90).max() > 1e-4 or not args.thickness:
            sys.exit("--cell-from needs an orthorhombic slab and --thickness")
        box = np.array([slab.cell[0, 0], slab.cell[1, 1], args.thickness])
    else:
        box = np.array(args.box)

    volume = np.prod(box)
    os.makedirs("packmol_files", exist_ok=True)
    blocks, solute_mass, excluded = [], 0.0, 0.0
    if args.solute:
        solute = read(args.solute)
        excluded = np.prod(np.ptp(solute.positions, axis=0) + 1.5)
        solute.positions += box / 2 - solute.positions.mean(0)
        write("packmol_files/solute.xyz", solute)
        blocks.append("structure packmol_files/solute.xyz\n  number 1\n"
                      "  fixed 0. 0. 0. 0. 0. 0.\nend structure")

    t = args.tolerance / 2
    inside = f"inside box {t:.3f} {t:.3f} {t:.3f} " + " ".join(f"{b - t:.3f}" for b in box)
    report = []
    other_mass = 0.0
    for k, spec in enumerate(args.molecules):
        name, mol = load_molecule(spec)
        mass = mol.get_masses().sum()
        if k == 0:
            continue
        n = int(spec.split(":")[1])
        other_mass += n * mass
        write(f"packmol_files/{name}.xyz", mol)
        blocks.append(f"structure packmol_files/{name}.xyz\n  number {n}\n  {inside}\n"
                      "end structure")
        report.append(f"{n} {name}")
    name0, mol0 = load_molecule(args.molecules[0])
    mass0 = mol0.get_masses().sum()
    free = volume - excluded
    n0 = int((args.density * free / (AMU_TO_G * 1e24) - other_mass) / mass0)
    if n0 <= 0:
        sys.exit("no room for the first molecule: check box, density and the added molecules")
    write(f"packmol_files/{name0}.xyz", mol0)
    blocks.insert(1 if args.solute else 0,
                  f"structure packmol_files/{name0}.xyz\n  number {n0}\n  {inside}\n"
                  "end structure")
    report.insert(0, f"{n0} {name0}")

    with open("packmol.inp", "w") as fh:
        seed = "" if args.seed is None else f"seed {args.seed}\n"
        fh.write(f"tolerance {args.tolerance}\nfiletype xyz\noutput {args.output}.xyz\n"
                 f"{seed}\n" + "\n\n".join(blocks) + "\n")
    print(f"box {np.round(box, 3)} A, free volume {free:.0f} A^3: " + ", ".join(report))
    print("-> packmol.inp")

    exe = shutil.which(args.packmol)
    if exe is None:
        print(f"packmol not found ('{args.packmol}'): run  packmol < packmol.inp  yourself")
        return
    with open("packmol.inp") as fin, open("packmol.log", "w") as flog:
        subprocess.run([exe], stdin=fin, stdout=flog, check=False)
    log = open("packmol.log").read()
    if "Success!" not in log:
        sys.exit("packmol did not converge (see packmol.log): lower the density or "
                 "enlarge the box")
    atoms = read(f"{args.output}.xyz")
    atoms.set_cell(box)
    atoms.pbc = True
    write(f"{args.output}.extxyz", atoms)
    print(f"packmol: Success! {len(atoms)} atoms -> {args.output}.extxyz")


if __name__ == "__main__":
    main()
