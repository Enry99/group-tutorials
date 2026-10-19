"""Small helpers shared by the scripts in tools/ (not meant to be run)."""
import io
import os
import re
import shlex
import sys

import numpy as np
from ase.io import read, write
from ase.io.espresso import read_fortran_namelist

from qe2xyz import read_pwi, strip

RY_TO_EV = 13.605693122994


def load_structure(path, index=-1):
    """Atoms from .pwi (any ibrav), .pwo (last frame by default) or any ASE format.

    A .pwo prints the cell of an scf/relax with 6 digits only (a(i) in units of
    alat). That is enough to look at it, not to build new inputs: a hexagonal
    cell becomes slightly non-hexagonal, and pw.x can then stop with
    "Error in routine d_matrix". If the .pwi with the same name exists and has
    the same cell to 1e-4 A (scf, relax), its exact cell is used.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in (".pwi", ".in"):
        return read_pwi(path)[0]
    if ext in (".pwo", ".out"):
        atoms = read(path, format="espresso-out", index=index)
        pwi = os.path.splitext(path)[0] + ".pwi"
        if os.path.exists(pwi):
            cell = read_pwi(pwi)[0].cell
            for frame in atoms if isinstance(atoms, list) else [atoms]:
                if np.allclose(frame.cell, cell, atol=1e-4):
                    frame.set_cell(cell)            # positions (Cartesian) unchanged
        return atoms
    return read(path, index=index)


def sort_z(atoms):
    """Group convention: atoms sorted along z (stable: ties keep their order)."""
    order = np.argsort(np.round(atoms.positions[:, 2], 4), kind="stable")
    return atoms[order]


def read_template(path):
    """Namelists and ATOMIC_SPECIES (label -> (mass, pseudo)) of a pw.x input."""
    with open(path) as fh:
        data, lines = read_fortran_namelist(fh)
    ntyp = int(data["system"]["ntyp"])
    species = {}
    for i, line in enumerate(lines):
        if strip(line).upper().startswith("ATOMIC_SPECIES"):
            for entry in lines[i + 1:]:
                s = strip(entry)
                if not s:
                    continue
                label, mass, pseudo = s.split()[:3]
                species[label] = (float(mass), pseudo)
                if len(species) == ntyp:
                    break
            break
    return data, species


def write_pwi(path, atoms, input_data, **kwargs):
    """ASE espresso-in writer + clean-up: a first comment line with the command
    that wrote the file, no empty namelists (&CELL, &FCP, &RISM), no '-0.0'."""
    buf = io.StringIO()
    write(buf, atoms, format="espresso-in", input_data=input_data, **kwargs)
    calc = input_data.get("control", {}).get("calculation", "scf")
    empty = "IONS|CELL|FCP|RISM" if calc in ("scf", "nscf", "bands") else "FCP|RISM"
    if calc == "relax":
        empty = "CELL|" + empty      # relax needs &IONS (even empty), not &CELL
    text = re.sub(rf"^&({empty})\n/\n", "", buf.getvalue(), flags=re.M)
    text = re.sub(r"(?<![\d.])-(0\.0+)(?!\d*[1-9])", r"\1", text)
    lines = [line.rstrip() for line in text.splitlines()]
    cmd = " ".join([os.path.basename(sys.argv[0])] + [shlex.quote(a) for a in sys.argv[1:]])
    with open(path, "w") as fh:
        fh.write(f"! written by: {cmd}\n" + "\n".join(lines) + "\n")


def kgrid_from_spacing(atoms, kspacing, slab=True):
    """Monkhorst-Pack grid for a k-spacing (1/A, 2pi included); n3 = 1 for slabs."""
    bnorm = 2.0 * np.pi * np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    grid = [max(1, int(np.ceil(b / kspacing))) for b in bnorm]
    if slab:
        grid[2] = 1
    return grid
