#!/usr/bin/env python3
"""Convert Quantum ESPRESSO files to extended XYZ (or any ASE format) for viewing.

Usage:
    python qe2xyz.py slab.pwi                 # -> slab.extxyz (any ibrav)
    python qe2xyz.py relax.pwo                # -> relax.extxyz, all ionic steps
    python qe2xyz.py relax.pwo --last         # only the final structure
    python qe2xyz.py neb.pwi                  # NEB input: one frame per given image
    python qe2xyz.py *.pwi -f cif             # other formats: cif, xsf, vasp, lammps-data
    python qe2xyz.py in.pwi -o check.extxyz

Then open the .extxyz in OVITO (basic version), ASE gui, VESTA (with -f cif), ...

Inputs (.pwi/.in) are read by this script, not by ASE, so every ibrav works
(cell built with the same formulas as QE's latgen.f90, checked against pw.x).
For inputs the extxyz also contains, per atom:
    qe_type    species index as in ATOMIC_SPECIES (1, 2, ...): colour by it to
               see Fe1/Fe2 (AFM) or other labels; the list is in the header
    if_pos     the three 0/1 flags (1 = free, 0 = fixed); default 1 1 1
    fixed      1 if all three flags are 0
    start_mag  starting_magnetization of the atom's species
Outputs (.pwo/.out) are read with ASE (energies, forces, magnetic moments kept).
Other extensions (.xsf, .xyz, POSCAR, ...) are passed to ASE; neb.axsf (NEB path,
all images) is read as xsf and wrapped into the cell.
Needs ase, numpy.

Bash alias (in ~/.bashrc), so it works without activating the environment:
    alias qe2xyz='/leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/python /leonardo/pub/userexternal/epedrett/tools/qe/qe2xyz.py'
"""
import argparse
import os
import re
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.io.espresso import label_to_symbol, read_fortran_namelist

BOHR = 0.529177210903  # Angstrom, as in QE 7.3 (constants.f90)
CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS",
         "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES", "ATOMIC_VELOCITIES",
         "HUBBARD", "SOLVENTS", "ADDITIONAL_K_POINTS", "TOTAL_CHARGE",
         "BEGIN_POSITIONS", "END_POSITIONS", "FIRST_IMAGE", "INTERMEDIATE_IMAGE",
         "LAST_IMAGE", "END_ENGINE_INPUT", "END")


def strip(line):
    return re.split(r"[!#]", line, maxsplit=1)[0].strip()


def unit_of(header):
    match = re.search(r"(alat|bohr|angstrom|crystal_sg|crystal)", header.lower())
    return match.group(1) if match else None


def celldm_from_system(system):
    """celldm(1..6) as QE builds it from celldm(i) or from A, B, C, cosAB, ..."""
    ibrav = int(system.get("ibrav", 0))
    celldm = np.zeros(7)  # 1-based
    has_celldm = any(f"celldm({i})" in system for i in range(1, 7))
    has_abc = "a" in system
    if has_celldm and has_abc:
        sys.exit("do not specify both celldm and a,b,c (pw.x would stop too)")
    if has_celldm:
        for i in range(1, 7):
            celldm[i] = float(system.get(f"celldm({i})", 0.0))
    elif has_abc:  # abc2celldm in latgen.f90
        a = float(system["a"])
        b = float(system.get("b", 0.0))
        c = float(system.get("c", 0.0))
        cosab = float(system.get("cosab", 0.0))
        cosac = float(system.get("cosac", 0.0))
        cosbc = float(system.get("cosbc", 0.0))
        celldm[1] = a / BOHR
        celldm[2] = b / a
        celldm[3] = c / a
        if ibrav in (14, 0):
            celldm[4], celldm[5], celldm[6] = cosbc, cosac, cosab
        elif ibrav in (-12, -13):
            celldm[5] = cosac
        elif ibrav in (-5, 5, 12, 13):
            celldm[4] = cosab
    return ibrav, celldm


def latgen(ibrav, d):
    """Lattice vectors in bohr, same conventions as QE latgen.f90 (d = celldm, 1-based)."""
    a1, a2, a3 = np.zeros(3), np.zeros(3), np.zeros(3)
    sr2, sr3 = np.sqrt(2.0), np.sqrt(3.0)
    if ibrav == 1:
        a1[0] = a2[1] = a3[2] = d[1]
    elif ibrav == 2:
        t = d[1] / 2
        a1[:] = [-t, 0, t]; a2[:] = [0, t, t]; a3[:] = [-t, t, 0]
    elif abs(ibrav) == 3:
        t = d[1] / 2
        a1[:] = a2[:] = a3[:] = t
        if ibrav < 0:
            a1[0] = -t; a2[1] = -t; a3[2] = -t
        else:
            a2[0] = -t; a3[0] = -t; a3[1] = -t
    elif ibrav == 4:
        a1[0] = d[1]; a2[0] = -d[1] / 2; a2[1] = d[1] * sr3 / 2; a3[2] = d[1] * d[3]
    elif abs(ibrav) == 5:
        t1, t2 = np.sqrt(1 + 2 * d[4]), np.sqrt(1 - d[4])
        if ibrav == 5:
            a2[1] = sr2 * d[1] * t2 / sr3; a2[2] = d[1] * t1 / sr3
            a1[0] = d[1] * t2 / sr2; a1[1] = -a1[0] / sr3; a1[2] = a2[2]
            a3[0] = -a1[0]; a3[1] = a1[1]; a3[2] = a2[2]
        else:
            a1[0] = d[1] * (t1 - 2 * t2) / 3; a1[1] = d[1] * (t1 + t2) / 3; a1[2] = a1[1]
            a2[:] = [a1[2], a1[0], a1[1]]
            a3[:] = [a1[1], a1[2], a1[0]]
    elif ibrav == 6:
        a1[0] = a2[1] = d[1]; a3[2] = d[1] * d[3]
    elif ibrav == 7:
        a2[0] = d[1] / 2; a2[1] = a2[0]; a2[2] = d[3] * d[1] / 2
        a1[:] = [a2[0], -a2[0], a2[2]]
        a3[:] = [-a2[0], -a2[0], a2[2]]
    elif ibrav == 8:
        a1[0] = d[1]; a2[1] = d[1] * d[2]; a3[2] = d[1] * d[3]
    elif abs(ibrav) == 9:
        a1[0] = 0.5 * d[1]
        if ibrav == 9:
            a1[1] = a1[0] * d[2]; a2[0] = -a1[0]; a2[1] = a1[1]
        else:
            a1[1] = -a1[0] * d[2]; a2[0] = a1[0]; a2[1] = -a1[1]
        a3[2] = d[1] * d[3]
    elif ibrav == 91:
        a1[0] = d[1]; a2[1] = d[1] * d[2] * 0.5; a2[2] = -d[1] * d[3] * 0.5
        a3[1] = a2[1]; a3[2] = -a2[2]
    elif ibrav == 10:
        a2[0] = 0.5 * d[1]; a2[1] = a2[0] * d[2]
        a1[0] = a2[0]; a1[2] = a2[0] * d[3]
        a3[1] = a2[0] * d[2]; a3[2] = a1[2]
    elif ibrav == 11:
        a1[:] = [0.5 * d[1], 0.5 * d[1] * d[2], 0.5 * d[1] * d[3]]
        a2[:] = [-a1[0], a1[1], a1[2]]
        a3[:] = [-a1[0], -a1[1], a1[2]]
    elif ibrav == 12:
        sen = np.sqrt(1 - d[4] ** 2)
        a1[0] = d[1]; a2[0] = d[1] * d[2] * d[4]; a2[1] = d[1] * d[2] * sen; a3[2] = d[1] * d[3]
    elif ibrav == -12:
        sen = np.sqrt(1 - d[5] ** 2)
        a1[0] = d[1]; a2[1] = d[1] * d[2]; a3[0] = d[1] * d[3] * d[5]; a3[2] = d[1] * d[3] * sen
    elif ibrav == 13:
        sen = np.sqrt(1 - d[4] ** 2)
        a1[0] = 0.5 * d[1]; a1[2] = -a1[0] * d[3]
        a2[0] = d[1] * d[2] * d[4]; a2[1] = d[1] * d[2] * sen
        a3[0] = a1[0]; a3[2] = -a1[2]
    elif ibrav == -13:
        sen = np.sqrt(1 - d[5] ** 2)
        a1[0] = 0.5 * d[1]; a1[1] = a1[0] * d[2]
        a2[0] = -a1[0]; a2[1] = a1[1]
        a3[0] = d[1] * d[3] * d[5]; a3[2] = d[1] * d[3] * sen
    elif ibrav == 14:
        singam = np.sqrt(1 - d[6] ** 2)
        term = 1 + 2 * d[4] * d[5] * d[6] - d[4] ** 2 - d[5] ** 2 - d[6] ** 2
        if term < 0:
            sys.exit("celldm do not make sense (ibrav = 14)")
        term = np.sqrt(term / (1 - d[6] ** 2))
        a1[0] = d[1]
        a2[0] = d[1] * d[2] * d[6]; a2[1] = d[1] * d[2] * singam
        a3[0] = d[1] * d[3] * d[5]
        a3[1] = d[1] * d[3] * (d[4] - d[5] * d[6]) / singam
        a3[2] = d[1] * d[3] * term
    else:
        sys.exit(f"ibrav = {ibrav} not supported")
    return np.array([a1, a2, a3])


def find_card(lines, name, start=0):
    pattern = re.compile(rf"^{name}(?![A-Z_])", re.I)   # ATOMIC_POSITIONS, ..._POSITIONS{crystal}
    for i in range(start, len(lines)):
        if pattern.match(strip(lines[i])):
            return i
    return None


def card_body(lines, i, n=None):
    """Data lines of the card starting at index i (comments removed)."""
    out = []
    for line in lines[i + 1:]:
        s = strip(line)
        if not s:
            continue
        if re.split(r"[\s({]", s, maxsplit=1)[0].upper() in CARDS:
            break
        out.append(s)
        if n is not None and len(out) == n:
            break
    return out


def read_pwi(path):
    with open(path) as fh:
        data, lines = read_fortran_namelist(fh)
    system = data["system"]
    nat, ntyp = int(system["nat"]), int(system["ntyp"])
    ibrav, celldm = celldm_from_system(system)

    # cell (Angstrom) and alat (bohr)
    if ibrav != 0:
        cell_bohr = latgen(ibrav, celldm)
        alat = celldm[1]
    else:
        i = find_card(lines, "CELL_PARAMETERS")
        if i is None:
            sys.exit("ibrav = 0 needs CELL_PARAMETERS")
        unit = unit_of(lines[i]) or ("alat" if celldm[1] > 0 else "bohr")
        vec = np.array([[float(x) for x in l.split()[:3]] for l in card_body(lines, i, 3)])
        if unit == "alat":
            if celldm[1] <= 0:
                sys.exit("CELL_PARAMETERS alat needs celldm(1) or A")
            cell_bohr = vec * celldm[1]
        elif unit == "angstrom":
            cell_bohr = vec / BOHR
        else:
            cell_bohr = vec
        alat = celldm[1] if celldm[1] > 0 else np.linalg.norm(cell_bohr[0])
    cell = cell_bohr * BOHR

    # species
    i = find_card(lines, "ATOMIC_SPECIES")
    labels = [l.split()[0] for l in card_body(lines, i, ntyp)]
    mags = [float(system.get(f"starting_magnetization({k + 1})", 0.0)) for k in range(ntyp)]

    # one ATOMIC_POSITIONS card (pw.x) or one per NEB image
    frames, start = [], 0
    while True:
        i = find_card(lines, "ATOMIC_POSITIONS", start)
        if i is None:
            break
        unit = unit_of(lines[i]) or "alat"
        if unit == "crystal_sg":
            sys.exit("ATOMIC_POSITIONS crystal_sg is not supported")
        body = card_body(lines, i, nat)
        if len(body) != nat:
            sys.exit(f"{path}: {len(body)} position lines, nat = {nat}")
        lab = [l.split()[0] for l in body]
        xyz = np.array([[float(x) for x in l.split()[1:4]] for l in body])
        flags = np.array([[int(x) for x in l.split()[4:7]] if len(l.split()) >= 7 else [1, 1, 1]
                          for l in body])
        if unit == "crystal":
            pos = xyz @ cell
        elif unit == "alat":
            pos = xyz * alat * BOHR
        elif unit == "bohr":
            pos = xyz * BOHR
        else:
            pos = xyz
        unknown = set(lab) - set(labels)
        if unknown:
            sys.exit(f"labels {sorted(unknown)} not in ATOMIC_SPECIES")
        atoms = Atoms([label_to_symbol(x) for x in lab], positions=pos, cell=cell, pbc=True)
        qe_type = np.array([labels.index(x) + 1 for x in lab])
        atoms.set_array("qe_type", qe_type)
        atoms.set_array("if_pos", flags)
        atoms.set_array("fixed", (flags.sum(axis=1) == 0).astype(int))
        atoms.set_array("start_mag", np.array([mags[t - 1] for t in qe_type]))
        atoms.info["qe_types"] = ",".join(labels)
        atoms.info["ibrav"] = ibrav
        frames.append(atoms)
        start = i + 1
    if not frames:
        sys.exit(f"{path}: no ATOMIC_POSITIONS card")
    return frames


def convert(path, args):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".pwi", ".in") or (ext not in (".pwo", ".out") and "&system" in
                                   open(path, errors="replace").read().lower()):
        frames = read_pwi(path)
        kind = "input"
    elif ext in (".pwo", ".out"):
        frames = read(path, format="espresso-out", index=-1 if args.last else ":",
                      results_required=False)
        frames = frames if isinstance(frames, list) else [frames]
        kind = "output"
    else:
        # neb.axsf of neb.x: ASE knows the format as "xsf" (not by extension, and
        # QE's leading spaces defeat the detection by content)
        fmt = "xsf" if ext == ".axsf" else None
        frames = read(path, format=fmt, index=-1 if args.last else ":")
        frames = frames if isinstance(frames, list) else [frames]
        if fmt == "xsf":
            for frame in frames:
                frame.wrap()     # neb.x writes minimum-image positions, some outside the cell
        kind = "ASE"
    out = args.output or os.path.splitext(path)[0] + "." + ("extxyz" if args.format == "extxyz"
                                                           else args.format.replace("-", "_"))
    if args.format == "extxyz":
        write(out, frames, format="extxyz")
    else:
        write(out, frames if len(frames) > 1 else frames[0], format=args.format)
    a = frames[-1]
    print(f"{path} ({kind}) -> {out}: {len(frames)} frame(s), {len(a)} atoms, "
          f"{a.get_chemical_formula()}, cell {np.round(a.cell.lengths(), 4)} A "
          f"{np.round(a.cell.angles(), 2)} deg")
    dmin = min_distance(a)
    if dmin < 0.7:
        print(f"  WARNING: shortest interatomic distance {dmin:.3f} A: overlapping atoms?")


def min_distance(atoms):
    if len(atoms) < 2:
        return np.inf
    from ase.neighborlist import neighbor_list
    d = neighbor_list("d", atoms, 3.0)
    return d.min() if len(d) else np.inf


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("-o", "--output", help="output name (only with one input file)")
    parser.add_argument("-f", "--format", default="extxyz", help="any ASE write format")
    parser.add_argument("--last", action="store_true", help="outputs: last frame only")
    args = parser.parse_args()
    if args.output and len(args.files) > 1:
        sys.exit("-o works with one input file")
    for path in args.files:
        convert(path, args)


if __name__ == "__main__":
    main()
