#!/usr/bin/env python3
"""Fit E(V) with the Birch-Murnaghan equation of state (ASE) and give a0.

Usage:
    python eos_fit.py scale_*.pwo
    python eos_fit.py scale_*.pwo --bravais cI     # also print the cubic a0

Reads volume and total energy from each .pwo, fits Birch-Murnaghan and prints
V0, E0, B0 and the isotropic scale factor s0 to apply to the template cell
(file names scale_<s>.pwo, as written by make_scan.py --scale).
--bravais cP | cF | cI converts V0 of the primitive cell to the conventional
cubic lattice parameter (a0^3 = 1, 4 or 2 primitive volumes).
Saves eos.png. Needs ase, numpy, matplotlib.
"""
import argparse
import re
import sys

import numpy as np
from ase.eos import EquationOfState
from ase.io import read
from ase.units import Bohr, kJ

FACTOR = {"cP": 1, "cF": 4, "cI": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--bravais", choices=FACTOR)
    parser.add_argument("--eos", default="birchmurnaghan",
                        help="any ASE EOS: birchmurnaghan, murnaghan, vinet, ...")
    parser.add_argument("--out", default="eos.png")
    args = parser.parse_args()

    vols, ens, refs = [], [], []
    for path in args.files:
        text = open(path, errors="replace").read()
        if "JOB DONE." not in text or "convergence NOT achieved" in text:
            print(f"warning: {path} not finished or SCF not converged, skipped")
            continue
        atoms = read(path, format="espresso-out", index=-1)
        vols.append(atoms.get_volume())
        ens.append(atoms.get_potential_energy())
        match = re.search(r"scale_(\d+\.\d+)\.pwo$", path)
        if match:
            refs.append(atoms.get_volume() / float(match.group(1)) ** 3)
    if len(vols) < 5:
        sys.exit("need at least 5 finished points for a meaningful fit")

    vols, ens = np.array(vols), np.array(ens)
    eos = EquationOfState(vols, ens, eos=args.eos)
    v0, e0, b = eos.fit()
    nat = len(atoms)
    print(f"EOS          : {args.eos}")
    print(f"V0           : {v0:.4f} A^3 per cell ({v0 / nat:.4f} A^3/atom)")
    print(f"E0           : {e0:.6f} eV per cell ({e0 / nat:.6f} eV/atom)")
    print(f"B0           : {b / kJ * 1.0e24:.1f} GPa")
    if refs:
        s0 = (v0 / np.mean(refs)) ** (1.0 / 3.0)
        print(f"scale s0     : {s0:.5f}  (multiply celldm(1)/A of the template by this)")
    if args.bravais:
        a0 = (v0 * FACTOR[args.bravais]) ** (1.0 / 3.0)
        print(f"a0 ({args.bravais})      : {a0:.4f} A  ({a0 / Bohr:.4f} bohr)")
    if not vols.min() < v0 < vols.max():
        print("WARNING: the minimum is outside the sampled volumes: shift the scan")

    import matplotlib
    matplotlib.use("Agg")
    eos.plot(args.out)
    print("saved", args.out)


if __name__ == "__main__":
    main()
