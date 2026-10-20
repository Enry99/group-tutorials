#!/usr/bin/env python3
"""Cut a hemispherical tip (flat face on top, curved apex down) from a bulk
crystal, optionally passivating the curved surface with H and/or OH.

Usage:
    python tip.py C_diamond.cif --radius 15
    python tip.py C_diamond.cif --radius 15 --passivate H:1.0 -o C_tip_H.extxyz
    python tip.py C_diamond.cif --radius 20 --passivate H:0.5 OH:0.5 --seed 3

What it does (based on create_clean.py, create_with_H_OH.py and pass.py by
E. Pedretti):
  1. repeats the bulk to contain a sphere of --radius A and keeps the atoms
     inside it and below its centre (lower hemisphere): the flat face at the
     top is where the load is applied, the apex touches the counter surface;
  2. removes atoms left with fewer than --min-cn neighbours (default 2:
     single dangling atoms are unphysical);
  3. --passivate X:f adds a passivant to a fraction f of the undercoordinated
     atoms of the curved surface (more than --flat-skip A below the flat
     face): one H or OH per atom, along the radial direction, at the sum of
     the covalent radii (O-H 0.97 A). Atoms missing two bonds get one
     passivant only; a passivant closer than --clash x (sum of covalent
     radii) to any other atom is skipped;
  4. puts the tip in a box with --vacuum A on each side (non periodic in
     practice) and sorts the atoms along z.
Neighbours: ASE natural cutoffs x 1.3 (as in pass.py). The flat-face area
pi R^2 is printed: use it as --load-area in interface.py, so that --load is
the nominal pressure on the flat face. Needs ase, numpy.
"""
import argparse
import os
import sys

import numpy as np
from ase import Atom
from ase.data import covalent_radii
from ase.io import read
from ase.neighborlist import natural_cutoffs, neighbor_list

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import extension, sort_z, write_structure  # noqa: E402

D_OH = 0.97


def coordination(atoms, mult):
    i = neighbor_list("i", atoms, natural_cutoffs(atoms, mult=mult))
    return np.bincount(i, minlength=len(atoms))


def clashes(atoms, host, new, mult):
    """True if a new atom is closer than mult * (sum of covalent radii) to any
    atom except its host (and the other atoms of the same passivant)."""
    others = np.delete(np.arange(len(atoms)), host)
    for at in new:
        d = np.linalg.norm(atoms.positions[others] - at.position, axis=1)
        if np.any(d < mult * (covalent_radii[atoms.numbers[others]] + covalent_radii[at.number])):
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bulk")
    parser.add_argument("--radius", type=float, required=True, help="A")
    parser.add_argument("--min-cn", type=int, default=2)
    parser.add_argument("--passivate", nargs="+", metavar="X:f", help="H and/or OH fractions")
    parser.add_argument("--flat-skip", type=float, default=3.0,
                        help="do not passivate within this distance (A) of the flat face")
    parser.add_argument("--mult", type=float, default=1.3)
    parser.add_argument("--clash", type=float, default=1.4,
                        help="skip a passivant closer than clash x (r_cov sum) to other atoms")
    parser.add_argument("--vacuum", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("-f", "--format", default="extxyz")
    parser.add_argument("--specorder", nargs="+")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)

    bulk = read(args.bulk)
    bulk.pbc = True
    cn_bulk = coordination(bulk * (2, 2, 2), args.mult).max()
    rep = [int(np.ceil(2 * args.radius / h)) + 1 for h in
           1.0 / np.linalg.norm(bulk.cell.reciprocal(), axis=1)]
    atoms = bulk * rep
    centre = atoms.cell[:].sum(0) / 2
    r = np.linalg.norm(atoms.positions - centre, axis=1)
    atoms = atoms[(r <= args.radius) & (atoms.positions[:, 2] <= centre[2])]
    atoms.pbc = False

    removed = 0
    while True:
        cn = coordination(atoms, args.mult)
        bad = cn < args.min_cn
        if not bad.any():
            break
        removed += bad.sum()
        atoms = atoms[~bad]
    print(f"hemisphere R = {args.radius} A: {len(atoms)} atoms "
          f"(removed {removed} with CN < {args.min_cn}); bulk CN {cn_bulk}")

    if args.passivate:
        fracs = {s.split(":")[0]: float(s.split(":")[1]) for s in args.passivate}
        if set(fracs) - {"H", "OH"} or sum(fracs.values()) > 1 + 1e-9:
            sys.exit("--passivate: H and/or OH, fractions summing to <= 1")
        cn = coordination(atoms, args.mult)
        ztop = atoms.positions[:, 2].max()
        sites = np.where((cn < cn_bulk) & (atoms.positions[:, 2] < ztop - args.flat_skip))[0]
        rng.shuffle(sites)
        start = 0
        for species, frac in fracs.items():
            chosen = sites[start:start + int(round(frac * len(sites)))]
            start += len(chosen)
            placed = 0
            for i in chosen:
                host = atoms[i]
                u = host.position - centre
                u /= np.linalg.norm(u)
                if species == "H":
                    new = [Atom("H", host.position + u * (covalent_radii[host.number]
                                                          + covalent_radii[1]))]
                else:
                    o = host.position + u * (covalent_radii[host.number] + covalent_radii[8])
                    new = [Atom("O", o), Atom("H", o + u * D_OH)]
                if clashes(atoms, i, new, args.clash):
                    continue
                for at in new:
                    atoms.append(at)
                placed += 1
            print(f"passivated {placed} of {len(sites)} undercoordinated surface atoms "
                  f"with {species} ({len(chosen) - placed} skipped: too close to others)")

    atoms.center(vacuum=args.vacuum)
    atoms.pbc = True
    atoms = sort_z(atoms)
    print(f"{len(atoms)} atoms, {atoms.get_chemical_formula()}; box "
          f"{np.round(atoms.cell.lengths(), 3)} A; flat-face area pi R^2 = "
          f"{np.pi * args.radius ** 2:.1f} A^2")
    stem = os.path.splitext(os.path.basename(args.bulk))[0]
    ext = extension(args.format)
    out = args.output or f"{stem}_tip_R{args.radius:g}.{ext}"
    write_structure(out, atoms, args.format, args.specorder)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
