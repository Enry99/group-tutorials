#!/usr/bin/env python3
"""Coherent interface between two crystals with the Zur-McGill lattice
matching (pymatgen ZSLGenerator + CoherentInterfaceBuilder).

Usage:
    python zsl_interface.py TiO2_rutile.cif 1 1 0 Fe_bcc.cif 1 1 0 --list
    python zsl_interface.py TiO2_rutile.cif 1 1 0 Fe_bcc.cif 1 1 0 --select 3 \\
        --film-layers 3 --substrate-layers 4

Arguments: film bulk and Miller index, substrate bulk and Miller index. The
film is strained onto the substrate: put as substrate the stiffer crystal,
or the one whose lattice parameter you want to keep.

What it does (based on ZurAlgo.ipynb by E. Pedretti):
  1. takes the conventional cells of both bulks and finds the supercells of
     the two surface lattices that match within the tolerances:
     --length-tol (fraction, default 0.03) on the lengths, --angle-tol
     (fraction, 0.01) on the angle, --area-tol (fraction, 0.09) on the
     areas, up to --max-area A^2;
  2. --list prints the terminations and the matches sorted by area: cell,
     strain on the film (von Mises), number of atoms with the requested
     thicknesses. Small area and small strain conflict: choose;
  3. --select i writes match i: the interface (<film>_<substrate>_<i>.extxyz,
     c along z, atoms sorted along z) and the two strained slabs separately
     (..._film.extxyz, ..._substrate.extxyz), e.g. for interface.py with
     molecules in between.
Thicknesses (--film-layers, --substrate-layers) are in unit planes; --gap
and --vacuum in A. Needs ase, numpy, pymatgen.
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read, write
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import sort_z  # noqa: E402


def conventional(path):
    s = AseAtomsAdaptor.get_structure(read(path))
    return SpacegroupAnalyzer(s).get_conventional_standard_structure()


def to_atoms(structure):
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.wrap(pbc=[True, True, False])   # pymatgen leaves the film outside the cell in x, y
    atoms.pbc = True
    atoms.positions[:, 2] -= atoms.positions[:, 2].min()
    return sort_z(atoms)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("film")
    parser.add_argument("film_miller", nargs=3, type=int)
    parser.add_argument("substrate")
    parser.add_argument("substrate_miller", nargs=3, type=int)
    parser.add_argument("--max-area", type=float, default=400.0)
    parser.add_argument("--length-tol", type=float, default=0.03)
    parser.add_argument("--angle-tol", type=float, default=0.01)
    parser.add_argument("--area-tol", type=float, default=0.09)
    parser.add_argument("--termination", type=int, default=0, help="index from --list")
    parser.add_argument("--film-layers", type=int, default=3)
    parser.add_argument("--substrate-layers", type=int, default=3)
    parser.add_argument("--gap", type=float, default=2.5)
    parser.add_argument("--vacuum", type=float, default=15.0)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--select", type=int)
    args = parser.parse_args()
    if not args.list and args.select is None:
        sys.exit("give --list or --select i")

    zsl = ZSLGenerator(max_area_ratio_tol=args.area_tol, max_area=args.max_area,
                       max_length_tol=args.length_tol, max_angle_tol=args.angle_tol,
                       bidirectional=True)
    builder = CoherentInterfaceBuilder(substrate_structure=conventional(args.substrate),
                                       film_structure=conventional(args.film),
                                       film_miller=tuple(args.film_miller),
                                       substrate_miller=tuple(args.substrate_miller),
                                       zslgen=zsl)
    terms = builder.terminations
    if not terms:
        sys.exit("no match: raise --max-area or the tolerances")
    print("terminations (film, substrate):")
    for k, t in enumerate(terms):
        print(f"  {k}: {t}")
    term = terms[args.termination]
    interfaces = list(builder.get_interfaces(term, gap=args.gap, vacuum_over_film=args.vacuum,
                                             film_thickness=args.film_layers,
                                             substrate_thickness=args.substrate_layers))
    # the same match can appear several times (equivalent supercells): keep one
    seen, unique = set(), []
    for itf in interfaces:
        a, b = itf.lattice.abc[:2]
        key = (round(a, 2), round(b, 2), round(itf.lattice.gamma, 1),
               round(itf.interface_properties["von_mises_strain"], 4), len(itf))
        if key not in seen:
            seen.add(key)
            unique.append(itf)
    unique.sort(key=lambda i: (i.lattice.a * i.lattice.b * np.sin(np.radians(i.lattice.gamma)),
                               i.interface_properties["von_mises_strain"]))
    print(f"termination {args.termination}: {len(interfaces)} matches, {len(unique)} distinct")
    if args.list:
        print("    i   a (A)    b (A)  gamma   area (A^2)  strain (%)  atoms")
        for k, itf in enumerate(unique):
            lat = itf.lattice
            area = lat.a * lat.b * np.sin(np.radians(lat.gamma))
            print(f"{k:5d} {lat.a:8.3f} {lat.b:8.3f} {lat.gamma:6.1f} {area:11.1f} "
                  f"{itf.interface_properties['von_mises_strain'] * 100:10.2f} {len(itf):6d}")
        return

    itf = unique[args.select]
    stem = (f"{os.path.splitext(os.path.basename(args.film))[0]}_"
            f"{os.path.splitext(os.path.basename(args.substrate))[0]}_{args.select}")
    atoms = to_atoms(itf)
    write(f"{stem}.extxyz", atoms)
    write(f"{stem}_film.extxyz", to_atoms(itf.film))
    write(f"{stem}_substrate.extxyz", to_atoms(itf.substrate))
    lat = itf.lattice
    print(f"match {args.select}: {len(atoms)} atoms ({len(itf.film)} film, "
          f"{len(itf.substrate)} substrate), a = {lat.a:.4f} b = {lat.b:.4f} A, "
          f"gamma = {lat.gamma:.2f} deg, film strain (von Mises) "
          f"{itf.interface_properties['von_mises_strain'] * 100:.2f} %")
    print(f"-> {stem}.extxyz, {stem}_film.extxyz, {stem}_substrate.extxyz")


if __name__ == "__main__":
    main()
