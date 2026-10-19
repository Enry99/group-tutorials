#!/usr/bin/env python3
"""Convert a structure file (POSCAR, CIF, ...) into a pw.x SCF input.

Typical use: structure downloaded from Materials Project.
    python struct_to_qe.py Fe2O3_mp-XXXX.cif --primitive  --kspacing 0.20 -o Fe2O3.pwi
    python struct_to_qe.py POSCAR --conventional --kspacing 0.20 -o cell.pwi
    python struct_to_qe.py POSCAR --kgrid 6 6 6 -o cell.pwi    # grid given by hand

What it does:
  * prints space group, number of atoms, and whether the cell is primitive;
  * optionally converts to the standard primitive or conventional cell (spglib);
  * writes ibrav = 0 with CELL_PARAMETERS angstrom and ATOMIC_POSITIONS crystal,
    pseudopotential names from group_pseudos.json (edit unknown ones by hand),
    and default &SYSTEM/&ELECTRONS settings (edit them for your task).
The output is a starting point: open it, check every line, visualize it.
Needs ase, numpy, spglib.
"""
import argparse
import json
import os
import sys

import numpy as np
import spglib
from ase import Atoms
from ase.data import atomic_numbers
from ase.io import read

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe_common import write_pwi  # noqa: E402


def to_spglib(atoms):
    return (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)


def from_spglib(cell):
    lattice, positions, numbers = cell
    return Atoms(numbers=numbers, cell=lattice, scaled_positions=positions, pbc=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure")
    cell_group = parser.add_mutually_exclusive_group()
    cell_group.add_argument("--primitive", action="store_true")
    cell_group.add_argument("--conventional", action="store_true")
    k_group = parser.add_mutually_exclusive_group()
    k_group.add_argument("--kspacing", type=float, default=0.20,
                         help="1/Angstrom with 2pi (default 0.20, then converge!)")
    k_group.add_argument("--kgrid", nargs=3, type=int)
    parser.add_argument("--symprec", type=float, default=1e-3)
    parser.add_argument("-o", "--output", default="structure.pwi")
    args = parser.parse_args()

    atoms = read(args.structure)
    dataset = spglib.get_symmetry_dataset(to_spglib(atoms), symprec=args.symprec)
    prim = spglib.standardize_cell(to_spglib(atoms), to_primitive=True,
                                   symprec=args.symprec)
    n_prim = len(prim[2])
    print(f"input            : {args.structure}, {len(atoms)} atoms, "
          f"{atoms.get_chemical_formula()}")
    print(f"space group      : {dataset.international} ({dataset.number})")
    print(f"primitive cell   : {n_prim} atoms -> input cell is "
          f"{'primitive' if n_prim == len(atoms) else f'{len(atoms) // n_prim}x primitive'}")

    if args.primitive:
        atoms = from_spglib(prim)
    elif args.conventional:
        atoms = from_spglib(spglib.standardize_cell(to_spglib(atoms), to_primitive=False,
                                                    symprec=args.symprec))
    print(f"written cell     : {len(atoms)} atoms, lengths {np.round(atoms.cell.lengths(), 4)} A, "
          f"angles {np.round(atoms.cell.angles(), 2)}")

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "group_pseudos.json")) as fh:
        group = json.load(fh)
    species = sorted(set(atoms.get_chemical_symbols()), key=lambda s: atomic_numbers[s])
    pseudos = {}
    for sym in species:
        pseudos[sym] = group.get(sym, f"{sym}.CHOOSE_ME.UPF")
        if sym not in group:
            print(f"WARNING: no group pseudo for {sym}: choose one (see 03_pseudopotentials.md)")

    if args.kgrid:
        kpts = tuple(args.kgrid)
    else:
        bnorm = 2.0 * np.pi * np.linalg.norm(atoms.cell.reciprocal(), axis=1)
        kpts = tuple(int(n) for n in np.maximum(1, np.ceil(bnorm / args.kspacing)))
    print(f"k-point grid     : {kpts}")

    input_data = {
        "control": {"calculation": "scf", "prefix": "pwscf", "outdir": "./WORK",
                    "pseudo_dir": "/leonardo/pub/userexternal/epedrett/pseudo",
                    "tprnfor": True, "tstress": True},
        "system": {"ecutwfc": 60, "ecutrho": 480, "occupations": "smearing",
                   "smearing": "gaussian", "degauss": 0.00735, "vdw_corr": "dft-d3",
                   "nspin": 1},
        "electrons": {"conv_thr": 1.0e-8, "mixing_beta": 0.3},
    }
    write_pwi(args.output, atoms, input_data,
              pseudopotentials=pseudos, kpts=kpts, crystal_coordinates=True)
    print(f"wrote {args.output}. Magnetic system? set nspin = 2 and starting_magnetization.")


if __name__ == "__main__":
    main()
