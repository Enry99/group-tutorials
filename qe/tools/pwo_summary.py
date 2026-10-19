#!/usr/bin/env python3
"""Quick health check of pw.x output files.

Usage:
    python pwo_summary.py run.pwo [more.pwo ...]

For each file prints: status, last total energy, SCF iterations of the last
step, total force, pressure, total/absolute magnetization, wall time.
Only the standard library is needed.

Status values:
    DONE       "JOB DONE." found and nothing went wrong
    SCF-FAIL   "convergence NOT achieved" (energies/forces are NOT usable)
    MAXSTEP    relax/md stopped at nstep before convergence
    TIMEOUT    max_seconds reached ("Maximum CPU time exceeded"): restart it
    CRASHED    no "JOB DONE." (killed by Slurm, error, or still running):
               look at the CRASH file and the Slurm .err file
A file with several runs (restart appended with ">>") gets the status of the last
run; energy, forces and pressure are always the last ones printed.
"""
import re
import sys

RY_TO_EV = 13.605693122994

PATTERNS = {
    "energy": re.compile(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", re.M),
    "final": re.compile(r"^\s*Final (?:energy|enthalpy)\s+=\s+(-?\d+\.\d+)\s+Ry", re.M),
    "niter": re.compile(r"convergence has been achieved in\s+(\d+)\s+iterations"),
    "force": re.compile(r"Total force =\s+(\d+\.\d+)"),
    "press": re.compile(r"total\s+stress.*P=\s*(-?\d+\.\d+)"),
    "totmag": re.compile(r"total magnetization\s+=\s+(-?\d+\.\d+)"),
    "absmag": re.compile(r"absolute magnetization\s+=\s+(-?\d+\.\d+)"),
    "nat": re.compile(r"number of atoms/cell\s+=\s+(\d+)"),
    "wall": re.compile(r"PWSCF\s+:.*CPU\s+(.*?)\s+WALL"),
}


def last(pattern, text):
    found = pattern.findall(text)
    return found[-1] if found else None


def summarize(path):
    with open(path, errors="replace") as fh:
        text = fh.read()

    # a restart appended with ">>" holds several runs: the status is that of the last one
    run = text.split("Program PWSCF v.")[-1]
    if "convergence NOT achieved" in run:
        status = "SCF-FAIL"
    elif "Maximum CPU time exceeded" in run:
        status = "TIMEOUT"
    elif "The maximum number of steps has been reached" in run:
        status = "MAXSTEP"
    elif "JOB DONE." in run:
        status = "DONE"
    else:
        status = "CRASHED"

    energy = last(PATTERNS["final"], text) or last(PATTERNS["energy"], text)
    nat = last(PATTERNS["nat"], text)
    out = {"file": path, "status": status}
    if energy is not None:
        e_ry = float(energy)
        out["E (Ry)"] = f"{e_ry:.8f}"
        if nat:
            out["E/atom (eV)"] = f"{e_ry * RY_TO_EV / int(nat):.5f}"
    for key, label in [("niter", "SCF it"), ("force", "F_tot (Ry/bohr)"),
                       ("press", "P (kbar)"), ("totmag", "M_tot (muB)"),
                       ("absmag", "M_abs (muB)"), ("wall", "wall")]:
        value = last(PATTERNS[key], text)
        if value is not None:
            out[label] = value.strip()
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for path in sys.argv[1:]:
        info = summarize(path)
        print(f"{info.pop('file')}: {info.pop('status')}")
        for key, value in info.items():
            print(f"    {key:16s} {value}")


if __name__ == "__main__":
    main()
