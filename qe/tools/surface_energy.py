#!/usr/bin/env python3
"""Surface energy of a slab from pw.x outputs.

Usage (example 06):
    python surface_energy.py Al111_6L.pwo Al_fcc_bulk.pwo

    Esurf = (E_slab - N * e_bulk) / (2 A)

E_slab: final energy of the relaxed slab (last energy in the .pwo).
e_bulk: energy per atom of the bulk, computed with the SAME pseudopotential,
        cutoffs, smearing, D3 and an equivalent k-point spacing.
A:      area of the surface cell. The factor 2 counts the two surfaces.
Valid when the two surfaces are equivalent (symmetric slab, all atoms relaxed).

Needs ase, numpy.
"""
import sys

import numpy as np
from ase.io import read

EV_A2_TO_J_M2 = 16.0217663


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    slab = read(sys.argv[1], format="espresso-out", index=-1)
    bulk = read(sys.argv[2], format="espresso-out", index=-1)
    e_slab = slab.get_potential_energy()
    e_bulk = bulk.get_potential_energy() / len(bulk)
    area = np.linalg.norm(np.cross(slab.cell[0], slab.cell[1]))
    gamma = (e_slab - len(slab) * e_bulk) / (2.0 * area)
    print(f"E_slab = {e_slab:.6f} eV ({len(slab)} atoms), e_bulk = {e_bulk:.6f} eV/atom")
    print(f"A = {area:.4f} A^2")
    print(f"gamma = {gamma:.4f} eV/A^2 = {gamma * EV_A2_TO_J_M2:.3f} J/m^2")


if __name__ == "__main__":
    main()
