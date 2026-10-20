"""Small helpers shared by the scripts in tools/ (not meant to be run)."""
import numpy as np
from ase.io import read, write

AMU_TO_G = 1.66053906660e-24      # g
A3_TO_CM3 = 1.0e-24               # cm^3
EV_PER_A_TO_N = 1.602176634e-9    # 1 eV/A in newton


def sort_z(atoms):
    """Group convention: atoms sorted along z (stable: ties keep their order)."""
    order = np.argsort(np.round(atoms.positions[:, 2], 4), kind="stable")
    return atoms[order]


def density(atoms):
    """Mass density in g/cm^3."""
    return atoms.get_masses().sum() * AMU_TO_G / (atoms.cell.volume * A3_TO_CM3)


def force_per_atom(pressure_gpa, area_a2, natoms):
    """Force (eV/A) on each of natoms atoms for a pressure (GPa) on an area (A^2)."""
    total_n = pressure_gpa * 1e9 * area_a2 * 1e-20
    return total_n / natoms / EV_PER_A_TO_N


def read_structure(path):
    """ASE read; a LAMMPS data file (.data, atom_style atomic) too, the
    elements taken from its masses."""
    if path.endswith(".data"):
        return read(path, format="lammps-data", atom_style="atomic")
    return read(path)


def extension(fmt):
    """File extension for an ASE format name (lammps-data -> data)."""
    return {"lammps-data": "data"}.get(fmt, fmt)


def write_structure(path, atoms, fmt="extxyz", specorder=None):
    """ASE write; lammps-data with masses, atom_style atomic and the type order
    of pair_coeff (specorder)."""
    if fmt == "lammps-data":
        kwargs = {"masses": True, "atom_style": "atomic"}
        if specorder:
            missing = set(atoms.get_chemical_symbols()) - set(specorder)
            if missing:
                raise SystemExit(f"--specorder misses {sorted(missing)}")
            kwargs["specorder"] = specorder
        write(path, atoms, format=fmt, **kwargs)
    else:
        write(path, atoms, format=fmt)
