#!/usr/bin/env python3
"""Random substitutional doping and interstitial atoms, e.g. for Si- or
H-doped amorphous carbon or hydrogenated silicon nitride.

Usage:
    python dope.py C_diamond.cif --repeat 5 4 4 --substitute Si:0.10 H:0.20
    python dope.py Si3N4.cif --repeat 2 2 4 --interstitial H:0.10 --seed 7
    python dope.py slab.extxyz --substitute N:0.05 --host C

What it does (based on dope_DLC.py and si3n4_doping.py by E. Pedretti):
  --substitute X:f  replaces a fraction f of ALL atoms with X, choosing among
                    the host atoms (--host, default: every species that is
                    not a dopant); several dopants are placed in turn;
  --interstitial X:f  adds f * N(start) atoms X at random positions, rejected
                    if closer than mult * (r_cov(X) + r_cov(neighbour)) to any
                    atom (--mult, default 1.4, covalent radii of ASE);
                    stops after --max-attempts random positions;
  --seed            makes the choice reproducible.
The result keeps the cell: set the density afterwards with density.py (H
atoms add little mass but a lot of volume in the real material). Needs ase,
numpy.
"""
import argparse
import os
import sys

import numpy as np
from ase import Atom
from ase.data import atomic_numbers, covalent_radii
from ase.geometry import get_distances
from ase.io import read

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import density, write_structure  # noqa: E402


def parse(items):
    out = {}
    for item in items or []:
        symbol, frac = item.split(":")
        out[symbol] = float(frac)
    return out


def insert(atoms, symbol, n, mult, rng, max_attempts):
    r_new = covalent_radii[atomic_numbers[symbol]]
    added, attempts = 0, 0
    while added < n and attempts < max_attempts:
        attempts += 1
        pos = rng.random(3) @ atoms.cell[:]
        d = get_distances(pos, atoms.positions, cell=atoms.cell, pbc=atoms.pbc)[1][0]
        if np.all(d >= mult * (covalent_radii[atoms.numbers] + r_new)):
            atoms.append(Atom(symbol, pos))
            added += 1
    if added < n:
        sys.exit(f"only {added} of {n} {symbol} placed after {attempts} attempts: "
                 "lower --mult or the fraction, or enlarge the cell")
    return atoms, attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structure")
    parser.add_argument("--repeat", nargs=3, type=int, default=(1, 1, 1))
    parser.add_argument("--substitute", nargs="+", metavar="X:f")
    parser.add_argument("--host", nargs="+", help="species that can be replaced")
    parser.add_argument("--interstitial", nargs="+", metavar="X:f")
    parser.add_argument("--mult", type=float, default=1.4)
    parser.add_argument("--max-attempts", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-f", "--format", default="extxyz")
    parser.add_argument("--specorder", nargs="+",
                        help="lammps-data: element order of the atom types (as in pair_coeff)")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    atoms = read(args.structure) * tuple(args.repeat)
    atoms.pbc = True
    n0 = len(atoms)
    subs, inters = parse(args.substitute), parse(args.interstitial)
    print(f"start: {atoms.get_chemical_formula()} ({n0} atoms), {density(atoms):.3f} g/cm^3")

    hosts = set(args.host) if args.host else set(atoms.get_chemical_symbols()) - set(subs)
    for symbol, frac in subs.items():
        candidates = [a.index for a in atoms if a.symbol in hosts]
        n = int(round(frac * n0))
        if n > len(candidates):
            sys.exit(f"{symbol}: {n} substitutions asked, {len(candidates)} host atoms")
        chosen = rng.choice(candidates, n, replace=False)
        atoms.symbols[chosen] = symbol
        print(f"substituted {n} host atoms with {symbol}")
    for symbol, frac in inters.items():
        n = int(round(frac * n0))
        atoms, attempts = insert(atoms, symbol, n, args.mult, rng, args.max_attempts)
        print(f"inserted {n} {symbol} ({attempts} attempts)")

    print(f"final: {atoms.get_chemical_formula()} ({len(atoms)} atoms), "
          f"{density(atoms):.3f} g/cm^3")
    stem = os.path.splitext(os.path.basename(args.structure))[0]
    tag = "_".join(f"{s}{round(f * 100)}" for s, f in {**subs, **inters}.items())
    ext = {"extxyz": "extxyz", "lammps-data": "data"}.get(args.format, args.format)
    out = args.output or f"{stem}_{tag}.{ext}"
    write_structure(out, atoms, args.format, args.specorder)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
