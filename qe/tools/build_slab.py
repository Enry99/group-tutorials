#!/usr/bin/env python3
"""Build a slab from the relaxed bulk and write the pw.x relax input.

Usage (example 06):
    python build_slab.py Al fcc111 --a0 4.012 --layers 6 \\
        --template ../03_lattice_parameter/01_eos/Al_fcc.pwi --kspacing 0.135

Arguments
  element surface   e.g. Al fcc111; surfaces: fcc100 fcc110 fcc111 bcc100
                    bcc110 bcc111 hcp0001 (ASE builders)
  --a0              bulk lattice parameter (A) from 03, same settings as the template
  --c               c (A) for hcp0001
  --layers          atomic layers
  --template        the bulk input: cutoffs, smearing, D3, pseudopotential and
                    SCF settings are copied from it, so slab and bulk match
  --kspacing        in-plane k-spacing (1/A, 2pi included) converged for the bulk
                    (tools/kgrid.py on a bulk .pwo prints it); n3 = 1
Options
  --vacuum 15       total vacuum (A), half on each side
  --repeat 1 1      lateral supercell (e.g. 2 2 for adsorption)
  --fixed 0         bottom layers to keep fixed (default 0: relax everything)
  --ecut W R        override ecutwfc/ecutrho (e.g. when an adsorbate needs more)
  --max-seconds N   override max_seconds (keep it ~5% below the job walltime)
  --magmom M        initial moment per atom (muB) for magnetic metals
  -o FILE           output name (default <element><surface>_<layers>L[_NxM].pwi)

The atoms are sorted along z (group convention). Writes also slab_info.json
with the parameters, used by add_adsorbate.py. Needs ase, numpy.
"""
import argparse
import json
import os
import sys

import ase.build
import numpy as np
from ase.constraints import FixAtoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe_common import kgrid_from_spacing, read_template, sort_z, write_pwi  # noqa: E402

SURFACES = ("fcc100", "fcc110", "fcc111", "bcc100", "bcc110", "bcc111", "hcp0001")
DROP_SYSTEM = {"ibrav", "a", "b", "c", "cosab", "cosac", "cosbc", "nat", "ntyp"}


def build(element, surface, a0, c, layers, repeat, vacuum):
    func = getattr(ase.build, surface)
    kwargs = {"a": a0, "vacuum": vacuum / 2.0, "orthogonal": False}
    if surface == "hcp0001":
        kwargs["c"] = c
    if surface in ("fcc100", "fcc110", "bcc100"):
        kwargs.pop("orthogonal")      # these builders have only the rectangular cell
    return func(element, size=(repeat[0], repeat[1], layers), **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("element")
    parser.add_argument("surface", choices=SURFACES)
    parser.add_argument("--a0", type=float, required=True)
    parser.add_argument("--c", type=float)
    parser.add_argument("--layers", type=int, required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--kspacing", type=float, required=True)
    parser.add_argument("--vacuum", type=float, default=15.0)
    parser.add_argument("--repeat", nargs=2, type=int, default=(1, 1))
    parser.add_argument("--fixed", type=int, default=0)
    parser.add_argument("--ecut", nargs=2, type=float, metavar=("WFC", "RHO"))
    parser.add_argument("--max-seconds", type=int)
    parser.add_argument("--magmom", type=float)
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    if args.surface == "hcp0001" and not args.c:
        sys.exit("hcp0001 needs --c")

    data, species = read_template(args.template)
    label = next((l for l in species if l.rstrip("0123456789") == args.element), None)
    if label is None:
        sys.exit(f"{args.element} not in ATOMIC_SPECIES of {args.template}")
    pseudo = species[label][1]

    slab = build(args.element, args.surface, args.a0, args.c, args.layers,
                 args.repeat, args.vacuum)
    # layer tags from ASE: 1 = top ... layers = bottom
    slab = sort_z(slab)
    if args.fixed:
        bottom = [a.index for a in slab if a.tag > args.layers - args.fixed]
        slab.set_constraint(FixAtoms(indices=bottom))
    magnetic = int(data["system"].get("nspin", 1)) == 2
    if magnetic:
        if args.magmom is None:
            sys.exit("the template has nspin = 2: give --magmom (muB per atom)")
        slab.set_initial_magnetic_moments([args.magmom] * len(slab))

    kpts = kgrid_from_spacing(slab, args.kspacing, slab=True)
    area = float(np.linalg.norm(np.cross(slab.cell[0], slab.cell[1])))
    nx, ny = args.repeat
    tag = "" if (nx, ny) == (1, 1) else f"_{nx}x{ny}"
    out = args.output or f"{args.element}{args.surface[3:]}_{args.layers}L{tag}.pwi"

    control = dict(data["control"])
    control.update({"calculation": "relax", "prefix": os.path.splitext(os.path.basename(out))[0],
                    "nstep": 200, "forc_conv_thr": 1.0e-3, "etot_conv_thr": 1.0e-5,
                    "disk_io": "low", "tprnfor": True})
    control.pop("tstress", None)        # the stress of a slab with vacuum is meaningless
    if args.max_seconds:
        control["max_seconds"] = args.max_seconds
    system = {k: v for k, v in data["system"].items()
              if k not in DROP_SYSTEM and not k.startswith(("celldm", "starting_magnetization"))}
    if args.ecut:
        system["ecutwfc"], system["ecutrho"] = (int(e) if e.is_integer() else e for e in args.ecut)
    electrons = dict(data["electrons"])
    electrons["mixing_mode"] = "local-TF"   # slabs: metal + vacuum
    electrons["mixing_beta"] = min(float(electrons.get("mixing_beta", 0.7)), 0.3)
    electrons["conv_thr"] = float(f"{1.0e-8 * len(slab):.1e}")  # total energy: scale with size
    input_data = {"control": control, "system": system, "electrons": electrons,
                  "ions": {"ion_dynamics": "bfgs"}}
    write_pwi(out, slab, input_data,
              pseudopotentials={args.element: pseudo}, kpts=kpts)

    info = {"element": args.element, "surface": args.surface, "a0": args.a0, "c": args.c,
            "layers": args.layers, "vacuum": args.vacuum, "repeat": [nx, ny],
            "fixed": args.fixed, "kpts": kpts, "magmom": args.magmom, "area": area,
            "template": args.template, "input": out}
    with open("slab_info.json", "w") as fh:
        json.dump(info, fh, indent=1)

    print(f"wrote {out}: {len(slab)} atoms in {args.layers} layers, "
          f"{args.fixed} fixed layer(s), sorted along z")
    print(f"  cell: {np.round(slab.cell.lengths(), 4)} A, angles {np.round(slab.cell.angles(), 2)}")
    print(f"  slab thickness {np.ptp(slab.positions[:, 2]):.3f} A, vacuum {args.vacuum} A, "
          f"area {area:.4f} A^2")
    print(f"  k-points {kpts[0]} {kpts[1]} 1, cutoffs {system['ecutwfc']}/{system['ecutrho']} Ry")


if __name__ == "__main__":
    main()
