#!/usr/bin/env python3
"""Write a series of pw.x inputs from one template, changing one parameter.

Usage (run it in the folder where the inputs should go):
    python make_scan.py template.pwi --ecut 30 40 50 60 70 80 --dual 8
    python make_scan.py template.pwi --kpts 4 6 8 10 12 14 16
    python make_scan.py template.pwi --scale 0.97 1.03 7
    python make_scan.py template.pwi --param degauss 0.002 0.005 0.00735 0.01

--ecut   sets ecutwfc and ecutrho = dual * ecutwfc (dual: 8 for US/PAW,
         4 for norm-conserving; see 00_basics/03_pseudopotentials.md).
--kpts   replaces the grid after "K_POINTS automatic" with n n n (shifts kept).
         Add --kz 1 for slabs (no sampling along the vacuum direction).
--scale  smin smax n: scales the lattice isotropically (celldm(1) or A, or
         CELL_PARAMETERS) by n factors between smin and smax. Positions in
         crystal or alat units follow automatically; angstrom/bohr positions
         are scaled too.
--param  any namelist variable, e.g. degauss, mixing_beta.
--set    NAME=VALUE ...: also change these variables in every input, e.g. an E(V)
         scan at another cutoff:  --scale 0.97 1.03 7 --set ecutwfc=90 ecutrho=720

Output names: ecut_060.pwi, k_08.pwi, scale_0.990.pwi, degauss_0.00735.pwi
Needs numpy.
"""
import argparse
import re
import sys

import numpy as np


def _pattern(name):
    """'name = value' at the start of a line or after a comma
    (ecutwfc = 60, ecutrho = 480)."""
    return re.compile(rf"(^|,)(\s*{re.escape(name)}\s*=\s*)([^,!\s]+)", re.I)


def set_namelist_value(text, name, value):
    """Replace 'name = old' with 'name = value' (comments after '!' left
    alone); error if not present once."""
    pattern, count, lines = _pattern(name), 0, []
    for line in text.split("\n"):
        code, bang, comment = line.partition("!")
        code, n = pattern.subn(lambda m: f"{m.group(1)}{m.group(2)}{value}", code)
        count += n
        lines.append(code + bang + comment)
    if count != 1:
        sys.exit(f"'{name}' found {count} times in the template (need exactly 1)")
    return "\n".join(lines)


def get_namelist_value(text, name):
    for line in text.split("\n"):
        match = _pattern(name).search(line.partition("!")[0])
        if match:
            return match.group(3).strip()
    return None


def set_kpoints(text, n, nz=None):
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if re.match(r"\s*K_POINTS\s*[({]?\s*automatic", line, re.I):
            old = lines[i + 1].split()
            shifts = old[3:6] if len(old) >= 6 else ["0", "0", "0"]
            kz = n if nz is None else nz
            lines[i + 1] = f"{n} {n} {kz}   {' '.join(shifts)}"
            return "\n".join(lines) + "\n"
    sys.exit("no 'K_POINTS automatic' card in the template")


def scale_lattice(text, s):
    """Isotropic scaling of the lattice (and of angstrom/bohr positions)."""
    ibrav = get_namelist_value(text, "ibrav")
    scaled = False
    for name in ("celldm(1)", "A"):
        value = get_namelist_value(text, name)
        if value is not None:
            new = float(value.replace("d", "e").replace("D", "e")) * s
            text = set_namelist_value(text, name, f"{new:.8f}")
            scaled = True
            break

    lines = text.splitlines()
    out = []
    mode = None
    for line in lines:
        head = line.strip().upper()
        if head.startswith("CELL_PARAMETERS"):
            mode = "cell_alat" if ("ALAT" in head or scaled) else "cell"
            out.append(line)
            continue
        if head.startswith("ATOMIC_POSITIONS"):
            mode = "pos" if ("ANGSTROM" in head or "BOHR" in head) else None
            out.append(line)
            continue
        if head.startswith(("K_POINTS", "ATOMIC_SPECIES", "HUBBARD", "&", "/",
                            "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES")):
            mode = None
        tokens = line.split()
        if mode == "cell" and len(tokens) == 3:
            vec = np.array([float(t) for t in tokens]) * s
            line = "  " + "  ".join(f"{v:14.8f}" for v in vec)
        elif mode == "pos" and len(tokens) >= 4:
            xyz = np.array([float(t) for t in tokens[1:4]]) * s
            line = f"{tokens[0]:4s}" + "".join(f"{v:16.10f}" for v in xyz)
            if len(tokens) > 4:
                line += "   " + " ".join(tokens[4:])
        out.append(line)
    if not scaled and int(ibrav or 0) != 0:
        sys.exit("ibrav != 0 but no celldm(1) or A found to scale")
    return "\n".join(out) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("template")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ecut", nargs="+", type=float)
    group.add_argument("--kpts", nargs="+", type=int)
    group.add_argument("--scale", nargs=3, type=float, metavar=("SMIN", "SMAX", "N"))
    group.add_argument("--param", nargs="+", metavar=("NAME", "VALUE"))
    parser.add_argument("--dual", type=float, default=8.0)
    parser.add_argument("--kz", type=int, default=None)
    parser.add_argument("--set", nargs="+", default=[], metavar="NAME=VALUE")
    args = parser.parse_args()

    with open(args.template) as fh:
        template = fh.read()
    for item in args.set:
        if "=" not in item:
            sys.exit(f"--set needs NAME=VALUE, got '{item}'")
        name, value = item.split("=", 1)
        template = set_namelist_value(template, name.strip(), value.strip())

    jobs = []
    if args.ecut:
        for ecut in args.ecut:
            text = set_namelist_value(template, "ecutwfc", f"{ecut:g}")
            text = set_namelist_value(text, "ecutrho", f"{ecut * args.dual:g}")
            jobs.append((f"ecut_{ecut:03.0f}.pwi", text))
    elif args.kpts:
        for k in args.kpts:
            jobs.append((f"k_{k:02d}.pwi", set_kpoints(template, k, args.kz)))
    elif args.scale:
        smin, smax, n = args.scale
        for s in np.linspace(smin, smax, int(n)):
            jobs.append((f"scale_{s:.3f}.pwi", scale_lattice(template, s)))
    else:
        name, values = args.param[0], args.param[1:]
        for value in values:
            jobs.append((f"{name}_{value}.pwi", set_namelist_value(template, name, value)))

    for fname, text in jobs:
        with open(fname, "w") as fh:
            fh.write(text)
        print("wrote", fname)


if __name__ == "__main__":
    main()
