#!/usr/bin/env python3
"""Cut slabs of any crystal and Miller index, one file per termination.

Usage:
    python slab.py TiO2_rutile.cif --miller 1 1 0 --layers 4 --list
    python slab.py TiO2_rutile.cif --miller 1 1 0 --layers 4 --vacuum 15
    python slab.py C_diamond.cif --miller 1 1 1 --layers 4 --keep-planes 8
    python slab.py Fe_bcc.cif --miller 1 1 0 --layers 6 \\
        --supercell 4 0 0 0 3 0 0 0 1 -f lammps-data --specorder Fe C O H

What it does (based on create_slab_new.py by E. Pedretti):
  1. reads the bulk (cif, POSCAR, .pwo, extxyz, ...) and takes its
     conventional standard cell (pymatgen SpacegroupAnalyzer);
  2. builds all the inequivalent terminations with pymatgen SlabGenerator;
     --layers is the thickness in unit planes (d_hkl);
  3. --keep-planes M keeps only the lowest M atomic planes: pymatgen may add
     planes so that top and bottom terminations are compatible; e.g. diamond
     (111) has 2 atomic planes per unit plane, so --layers 4 --keep-planes 8;
  4. makes c perpendicular to the surface, sets the vacuum (A) along z and
     applies the in-plane supercell (rows = new vectors, ASE/pymatgen
     convention; see qe/06_slab_al111, "Orthorhombic cells");
  5. sorts the atoms along z (group convention) and writes one file per
     termination: <stem>_<hkl>_t<i>.extxyz (or another ASE format, -f).
--list prints the terminations (label, atoms, formula, stoichiometry, top and
bottom plane) without writing. Needs ase, numpy, pymatgen.
"""
import argparse
import os
import sys

import numpy as np
from ase.geometry import get_layers
from ase.neighborlist import neighbor_list
from ase.io import read
from pymatgen.core.interface import label_termination
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import extension, sort_z, write_structure  # noqa: E402


def planes(atoms, tol):
    """Index of the atomic plane (along z) of each atom, from the bottom."""
    return get_layers(atoms, (0, 0, 1), tolerance=tol)[0]


def min_distance(atoms, rmax=4.0):
    d = neighbor_list("d", atoms, rmax)
    return d.min() if len(d) else rmax


def unwrap_z(atoms):
    """Make the slab contiguous along z (pymatgen may wrap some atoms through
    the vacuum) and put its lowest atom at z = 0."""
    sz = atoms.get_scaled_positions(wrap=True)[:, 2]
    order = np.sort(sz)
    gaps = np.diff(np.append(order, order[0] + 1.0))
    start = order[(np.argmax(gaps) + 1) % len(order)]   # first atom after the vacuum
    frac = atoms.get_scaled_positions(wrap=True)
    frac[:, 2] = (frac[:, 2] - start) % 1.0
    atoms.set_scaled_positions(frac)
    atoms.positions[:, 2] -= atoms.positions[:, 2].min()
    return sort_z(atoms)


def plane_formula(atoms, layer, which):
    idx = np.where(layer == (layer.max() if which == "top" else 0))[0]
    return atoms[idx].get_chemical_formula(mode="hill")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bulk")
    parser.add_argument("--miller", nargs=3, type=int, metavar=("H", "K", "L"), required=True)
    parser.add_argument("--layers", type=int, required=True, help="thickness in unit planes")
    parser.add_argument("--keep-planes", type=int, metavar="M",
                        help="keep only the lowest M atomic planes")
    parser.add_argument("--vacuum", type=float, default=15.0, help="vacuum along z (A)")
    parser.add_argument("--supercell", nargs=9, type=int, metavar="P",
                        help="in-plane supercell matrix, 9 integers by rows")
    parser.add_argument("--symmetrize", action="store_true",
                        help="equivalent top and bottom surfaces (may remove atoms)")
    parser.add_argument("--ztol", type=float, default=0.1,
                        help="tolerance (A) to group atoms into planes")
    parser.add_argument("--list", action="store_true", help="print the terminations only")
    parser.add_argument("-f", "--format", default="extxyz")
    parser.add_argument("--specorder", nargs="+",
                        help="lammps-data: element order of the atom types (as in pair_coeff)")
    args = parser.parse_args()

    structure = AseAtomsAdaptor.get_structure(read(args.bulk))
    conventional = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    bulk_formula = conventional.composition.reduced_composition
    d_bulk = min_distance(AseAtomsAdaptor.get_atoms(conventional))

    gen = SlabGenerator(conventional, tuple(args.miller), min_slab_size=args.layers,
                        min_vacuum_size=3, in_unit_planes=True, center_slab=False,
                        primitive=True, max_normal_search=max(map(abs, args.miller)) + 2)
    slabs = gen.get_slabs(symmetrize=args.symmetrize)
    hkl = "".join(str(i) for i in args.miller)
    if not slabs:
        sys.exit("no slab generated (with --symmetrize, try without)")

    stem = os.path.splitext(os.path.basename(args.bulk))[0]
    print(f"bulk {conventional.composition.reduced_formula}, conventional cell "
          f"{np.round(conventional.lattice.abc, 4)} A; {len(slabs)} termination(s) for {hkl}")
    for i, slab in enumerate(slabs):
        label = label_termination(slab)
        slab = slab.get_orthogonal_c_slab()
        atoms = sort_z(AseAtomsAdaptor.get_atoms(slab))
        layer = planes(atoms, args.ztol)
        nplanes = layer.max() + 1
        if args.keep_planes:
            if args.keep_planes > nplanes:
                sys.exit(f"t{i}: only {nplanes} atomic planes, cannot keep {args.keep_planes}")
            atoms = atoms[layer < args.keep_planes]
            layer = layer[layer < args.keep_planes]

        # c along z, vacuum on top (atoms start at z = 0)
        atoms = unwrap_z(atoms)
        thickness = np.ptp(atoms.positions[:, 2])
        cell = atoms.cell[:]
        cell[2] = [0.0, 0.0, thickness + args.vacuum]
        atoms.set_cell(cell)
        if args.supercell:
            from ase.build import make_supercell
            atoms = sort_z(make_supercell(atoms, np.reshape(args.supercell, (3, 3))))
            layer = planes(atoms, args.ztol)
        atoms.wrap()

        comp = AseAtomsAdaptor.get_structure(atoms).composition
        stoich = "stoichiometric" if comp.reduced_composition == bulk_formula else \
            f"NOT stoichiometric ({comp.reduced_formula} vs {bulk_formula.reduced_formula})"
        d_min = min_distance(atoms)
        if d_min < 0.9 * d_bulk:
            print(f"    WARNING: shortest distance {d_min:.3f} A, bulk {d_bulk:.3f} A: "
                  "overlapping atoms, do not use this slab")
        a, b = atoms.cell.lengths()[:2]
        gamma = atoms.cell.angles()[2]
        print(f"t{i} {label}: {len(atoms)} atoms, {nplanes} atomic planes"
              f"{f' (kept {args.keep_planes})' if args.keep_planes else ''}, {stoich}")
        print(f"    bottom plane {plane_formula(atoms, layer, 'bottom')}, "
              f"top plane {plane_formula(atoms, layer, 'top')}; "
              f"a = {a:.4f} b = {b:.4f} A, gamma = {gamma:.2f} deg, "
              f"thickness {thickness:.3f} A")
        if not args.list:
            ext = extension(args.format)
            out = f"{stem}_{hkl}_t{i}.{ext}"
            write_structure(out, atoms, args.format, args.specorder)
            print(f"    -> {out}")


if __name__ == "__main__":
    main()
