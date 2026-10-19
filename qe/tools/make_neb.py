#!/usr/bin/env python3
"""Write a neb.x input from two RELAXED end points.

Usage (example 08, from 08_neb_h_al111/02_neb):
    python make_neb.py ../../07_adsorption_h_al111/03_slab_h/Al111_4L_2x2_H_fcc.pwo \\
                       ../01_second_fcc/Al111_4L_2x2_H_fcc_10.pwo \\
                       --via ../../07_adsorption_h_al111/03_slab_h/Al111_4L_2x2_H_hcp.pwo --images 6

Arguments
  initial_pwo, final_pwo   relaxed end points; the .pwi of the initial one (same
                           name) gives all engine settings, species and fixed atoms
Options
  --images 6        num_of_images INCLUDING the two end points (6 = 4 moving
                    images = one per GPU with jobscripts/neb.sh)
  --via PWO         a structure the initial path must pass through (INTERMEDIATE_IMAGE);
                    repeat for more. neb.x spreads the images along the broken line.
  --prefix neb      prefix of the NEB files (neb.path, neb.dat, neb.axsf, ...)
  --max-seconds N   default 85000 (24 h walltime)

The atoms must be in the same order in both end points. Never sort or
renumber one end point alone: build both from the same structure (as 07 does)
and keep the order pw.x wrote. The script stops if the element sequence differs
and warns if more than one atom moves by more than 1.5 A (minimum image).
Stage 1 is written: no climbing image, loose path_thr. Stage 2: see README 08.
Writes neb.pwi. Needs ase, numpy.
"""
import argparse
import os
import sys

import numpy as np
from ase.io.espresso import read_fortran_namelist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe2xyz import read_pwi, strip  # noqa: E402
from qe_common import load_structure  # noqa: E402

SKIP = {"calculation", "nstep", "forc_conv_thr", "etot_conv_thr"}


def fortran(value):
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return f"'{value}'"
    if isinstance(value, float):
        return f"{value:.6g}".replace("e", "d")
    return str(value)


def namelist(name, data, skip=(), comments=None):
    comments = comments or {}
    lines = [f"&{name.upper()}"]
    for k, v in data.items():
        if k in skip:
            continue
        line = f"  {k} = {fortran(v)}"
        if k in comments:
            line = f"{line:34s}! {comments[k]}"
        lines.append(line)
    return "\n".join(lines + ["/"])


PATH_COMMENTS = {
    "restart_mode": "stage 2: 'restart'",
    "num_of_images": "including the two end points",
    "CI_scheme": "stage 2: 'auto' (climbing image)",
    "path_thr": "eV/A; stage 2: 0.05",
    "k_max": "spring constants (a.u.), QE default",
}


def card(lines, header):
    """Header line + data lines of the first card called `header`."""
    for i, line in enumerate(lines):
        if strip(line).upper().startswith(header):
            out = [strip(line)]
            for entry in lines[i + 1:]:
                s = strip(entry)
                if not s:
                    continue
                if s.split()[0].upper() in ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS",
                                            "CELL_PARAMETERS", "HUBBARD"):
                    break
                out.append(s)
            return out
    sys.exit(f"card {header} not found in the template")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("initial_pwo")
    parser.add_argument("final_pwo")
    parser.add_argument("--images", type=int, default=6)
    parser.add_argument("--prefix", default="neb")
    parser.add_argument("--max-seconds", type=int, default=85000)
    parser.add_argument("--via", action="append", default=[], metavar="PWO")
    args = parser.parse_args()

    template = args.initial_pwo[:-4] + ".pwi"
    with open(template) as fh:
        data, lines = read_fortran_namelist(fh)
    ref = read_pwi(template)[0]
    labels = card(lines, "ATOMIC_POSITIONS")[1:]
    labels = [l.split()[0] for l in labels][:len(ref)]
    flags = ref.arrays["if_pos"]

    paths = [args.initial_pwo] + args.via + [args.final_pwo]
    images = []
    for path in paths:
        if "End of BFGS Geometry Optimization" not in open(path, errors="replace").read():
            sys.exit(f"{path}: relaxation not finished")
        atoms = load_structure(path)
        if list(atoms.get_chemical_symbols()) != list(ref.get_chemical_symbols()):
            sys.exit(f"{path}: element sequence differs from {template} (atom order!)")
        images.append(atoms)

    for a, b, pa, pb in zip(images[:-1], images[1:], paths[:-1], paths[1:]):
        d = b.positions - a.positions
        frac = np.linalg.solve(a.cell.T, d.T).T
        frac -= np.round(frac)                          # minimum image
        disp = np.linalg.norm(frac @ a.cell, axis=1)
        i = int(disp.argmax())
        print(f"{os.path.basename(pa)} -> {os.path.basename(pb)}: largest displacement "
              f"{disp[i]:.3f} A (atom {i + 1}, {labels[i]})")
        if np.sort(disp)[-2] > 1.5:
            print("WARNING: more than one atom moves by > 1.5 A: is the atom order the same?")

    cmd = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])
    out = [f"! written by: {cmd}", "! template with comments: 08_neb_h_al111/neb_template.pwi",
           "BEGIN", "BEGIN_PATH_INPUT"]
    out.append(namelist("path", {
        "string_method": "neb",
        "restart_mode": "from_scratch",   # stage 2: 'restart'
        "nstep_path": 200,
        "num_of_images": args.images,
        "opt_scheme": "broyden",
        "CI_scheme": "no-CI",             # stage 2: 'auto' (climbing image)
        "path_thr": 0.1,                  # eV/A; stage 2: 0.05
        "k_max": 0.1, "k_min": 0.1,       # spring constants (a.u.), QE default
        "ds": 1.0,
        "minimum_image": True,
    }, comments=PATH_COMMENTS))
    out += ["END_PATH_INPUT", "", "BEGIN_ENGINE_INPUT"]
    control = dict(data["control"], prefix=args.prefix, max_seconds=args.max_seconds)
    out.append(namelist("control", control, SKIP))
    out.append(namelist("system", data["system"]))
    out.append(namelist("electrons", data["electrons"]))
    out.append("")
    for header in ("ATOMIC_SPECIES", "K_POINTS", "CELL_PARAMETERS"):
        out += card(lines, header) + [""]
    out.append("BEGIN_POSITIONS")
    names = ["FIRST_IMAGE"] + ["INTERMEDIATE_IMAGE"] * len(args.via) + ["LAST_IMAGE"]
    for name, atoms in zip(names, images):
        out += [name, "ATOMIC_POSITIONS angstrom"]
        for k, (label, pos) in enumerate(zip(labels, atoms.positions)):
            fix = "" if flags[k].all() else "   " + " ".join(str(f) for f in flags[k])
            out.append(f"{label:4s} {pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}{fix}")
    out += ["END_POSITIONS", "END_ENGINE_INPUT", "END"]
    with open("neb.pwi", "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"wrote neb.pwi: {args.images} images ({args.images - 2} moving), {len(ref)} atoms, "
          f"prefix '{args.prefix}'")


if __name__ == "__main__":
    main()
