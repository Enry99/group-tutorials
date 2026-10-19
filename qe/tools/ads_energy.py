#!/usr/bin/env python3
"""Adsorption energy from three pw.x outputs.

    E_ads = E(slab+X) - E(slab) - n * E(X)

X is the adsorbate as an isolated species (atom or molecule in a box), n the
number of X units on the slab (default 1). Negative E_ads = bound.

Usage (example 07, from 07_adsorption_h_al111):
    python ads_energy.py 03_slab_h/Al111_4L_2x2_H_fcc.pwo 02_clean_slab/Al111_4L_2x2.pwo \\
                         01_gas_phase/H_atom.pwo
    # dissociative adsorption of H2 (two H on the slab, H2 as reference):
    python ads_energy.py 03_slab_h/Al111_4L_2x2_H_fcc_H_fcc_11.pwo 02_clean_slab/Al111_4L_2x2.pwo \\
                         01_gas_phase/H2_box.pwo

Checks: the three runs must share pseudopotentials, cutoffs, smearing and D3
(the script compares the cutoffs in the outputs), and the atoms must add up:
N(slab+X) - N(slab) = n * N(X).
Needs ase.
"""
import argparse
import re

from ase.io import read


def cutoffs(path):
    text = open(path, errors="replace").read()
    wfc = re.search(r"kinetic-energy cutoff\s+=\s+([\d.]+)", text)
    rho = re.search(r"charge density cutoff\s+=\s+([\d.]+)", text)
    return (wfc.group(1) if wfc else "?", rho.group(1) if rho else "?")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("system", help="slab + adsorbate (.pwo)")
    parser.add_argument("slab", help="clean slab (.pwo)")
    parser.add_argument("reference", help="isolated adsorbate (.pwo)")
    parser.add_argument("-n", "--nref", type=float, default=1.0,
                        help="number of adsorbate units on the slab (default 1)")
    args = parser.parse_args()

    energies, natoms = {}, {}
    for label, path in [("slab+X", args.system), ("slab", args.slab), ("X", args.reference)]:
        atoms = read(path, format="espresso-out", index=-1)
        energies[label], natoms[label] = atoms.get_potential_energy(), len(atoms)
        print(f"{label:7s} {energies[label]:16.6f} eV  ({len(atoms)} atoms, "
              f"ecutwfc/ecutrho = {'/'.join(cutoffs(path))} Ry)  {path}")
    if len({cutoffs(p) for p in (args.system, args.slab, args.reference)}) > 1:
        print("WARNING: different cutoffs: the energies cannot be combined")
    if natoms["slab+X"] - natoms["slab"] != args.nref * natoms["X"]:
        print(f"WARNING: atoms do not add up: {natoms['slab+X']} - {natoms['slab']} "
              f"!= {args.nref:g} x {natoms['X']} (wrong reference or --nref?)")
    e_ads = energies["slab+X"] - energies["slab"] - args.nref * energies["X"]
    print(f"\nE_ads = {e_ads:.4f} eV")


if __name__ == "__main__":
    main()
