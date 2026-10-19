#!/usr/bin/env python3
"""Put one atom on a RELAXED slab and write the slab+adsorbate relax input.

Usage (example 07, from 03_slab_h):
    python add_adsorbate.py ../02_clean_slab/Al111_4L_2x2.pwo H fcc --height 1.0

Arguments
  slab_pwo   output of the slab relaxation (clean, or already with adsorbates);
             its .pwi (same name) and the slab_info.json written by
             build_slab.py must be in the same folder
  element    adsorbate atom (H, O, ...)
  site       ASE site names: fcc111 / hcp0001: ontop bridge fcc hcp;
             bcc110: ontop longbridge shortbridge hollow; fcc100: ontop bridge hollow
Options
  --height   start height above the top layer (A); the relax finds the final one
  --pseudo   pseudopotential of the adsorbate (default: tools/group_pseudos.json)
  --magmom   initial moment of the adsorbate (muB), only for nspin = 2
  --offset I J  same site type, I and J surface cells further along a1 and a2
             (e.g. the neighbouring fcc site for an fcc -> fcc NEB, example 08)

All pw.x settings come from the slab input, so slab and slab+adsorbate are
computed identically (cutoffs, k-points, smearing, D3, constraints). The atoms
are sorted along z: the adsorbate, on top, comes last. Run it again on the
relaxed slab+adsorbate for a second adsorbate (example 07, H2 exercise).
Writes <slab prefix>_<element>_<site>.pwi (+ _<I><J> with --offset). Needs ase, numpy.
"""
import argparse
import json
import os
import sys

import numpy as np
from ase import Atom
from ase.constraints import FixAtoms
from ase.io.espresso import read_fortran_namelist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe2xyz import read_pwi  # noqa: E402
from qe_common import load_structure, read_template, sort_z, write_pwi  # noqa: E402
from build_slab import build  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("slab_pwo")
    parser.add_argument("element")
    parser.add_argument("site")
    parser.add_argument("--height", type=float, default=1.0)
    parser.add_argument("--pseudo")
    parser.add_argument("--magmom", type=float, default=0.0)
    parser.add_argument("--offset", nargs=2, type=int, default=(0, 0), metavar=("I", "J"))
    args = parser.parse_args()

    folder = os.path.dirname(os.path.abspath(args.slab_pwo))
    with open(os.path.join(folder, "slab_info.json")) as fh:
        info = json.load(fh)
    pwi = args.slab_pwo[:-4] + ".pwi"
    text = open(args.slab_pwo, errors="replace").read()
    if "End of BFGS Geometry Optimization" not in text:
        sys.exit(f"{args.slab_pwo}: the clean-slab relaxation did not finish")

    slab = load_structure(args.slab_pwo)              # relaxed positions, same order as .pwi
    reference = read_pwi(pwi)[0]                      # for the fixed-atom flags
    if len(slab) != len(reference):
        sys.exit(f"{args.slab_pwo} and {pwi} have different numbers of atoms")

    # site position from the ideal slab: ASE stores the sites in the 1x1 surface cell
    ideal = build(info["element"], info["surface"], info["a0"], info["c"],
                  info["layers"], info["repeat"], info["vacuum"])
    sites = ideal.info["adsorbate_info"]["sites"]
    if args.site not in sites:
        sys.exit(f"site '{args.site}' not available for {info['surface']}: {sorted(sites)}")
    cell2d = np.array(ideal.info["adsorbate_info"]["cell"])
    # ASE builders put a top-layer atom at the origin and give the sites relative
    # to it. Same reference in the relaxed slab: the top-layer atom nearest to the
    # origin (minimum image in the plane)
    # (top layer = slab element only, so the slab may already carry adsorbates)
    host = np.array(slab.get_chemical_symbols()) == info["element"]
    ztop = slab.positions[host, 2].max()
    top = np.where(host & (slab.positions[:, 2] > ztop - 0.5))[0]
    frac = slab.cell.scaled_positions(slab.positions[top])[:, :2]
    frac -= np.round(frac)
    d = np.dot(frac, slab.cell[:2, :2])               # in-plane vectors from the origin
    xy0 = d[np.argmin(np.linalg.norm(d, axis=1))]
    xy = xy0 + np.dot(np.array(sites[args.site]) + np.array(args.offset), cell2d)
    z = slab.positions[top, 2].mean() + args.height
    slab.set_constraint()
    slab.set_tags(0)
    slab.append(Atom(args.element, (xy[0], xy[1], z), tag=1))   # tag: find it after sorting
    slab.wrap()

    fixed = [i for i, f in enumerate(reference.arrays["fixed"]) if f]
    if fixed:
        slab.set_constraint(FixAtoms(indices=fixed))
    with open(pwi) as fh:
        input_data, _ = read_fortran_namelist(fh)
    _, species = read_template(pwi)
    pseudos = {l.rstrip("0123456789"): p for l, (m, p) in species.items()}
    if args.pseudo:
        pseudos[args.element] = args.pseudo
    elif args.element not in pseudos:
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "group_pseudos.json")) as fh:
            pseudos[args.element] = json.load(fh)[args.element]
    if int(input_data["system"].get("nspin", 1)) == 2:
        moments = [info["magmom"]] * (len(slab) - 1) + [args.magmom]
        slab.set_initial_magnetic_moments(moments)
        for key in list(input_data["system"]):
            if key.startswith("starting_magnetization"):
                del input_data["system"][key]
    slab = sort_z(slab)

    prefix = f"{input_data['control']['prefix']}_{args.element}_{args.site}"
    if tuple(args.offset) != (0, 0):
        prefix += f"_{args.offset[0]}{args.offset[1]}"
    input_data["control"]["prefix"] = prefix
    out = prefix + ".pwi"
    write_pwi(out, slab, input_data,
              pseudopotentials=pseudos, kpts=tuple(info["kpts"]))
    i_ads = int(np.where(slab.get_tags() == 1)[0][0])
    ads = slab[i_ads]
    others = [i for i in range(len(slab)) if i != i_ads]
    dmin = slab.get_distances(i_ads, others, mic=True).min()
    print(f"wrote {out}: {args.element} on {args.site} at z = {ads.position[2]:.3f} A, "
          f"nearest atom {dmin:.3f} A; {len(slab)} atoms sorted along z")


if __name__ == "__main__":
    main()
