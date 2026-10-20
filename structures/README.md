# Structure generation

Scripts and examples to build the starting structures of our simulations:
slabs, coherent interfaces, amorphous and doped bulks, sliding interfaces
with a lubricant, tips. They do not depend on the code: the results are
`.extxyz` files (for OVITO, QE, any ASE calculator) or LAMMPS data files.

## Contents

| # | Folder | You learn |
|---|---|---|
| 1 | `01_slabs` | slabs of any crystal and Miller index, terminations, orthorhombic cells |
| 2 | `02_coherent_interface` | heterointerfaces with lattice matching (Zur-McGill) |
| 3 | `03_amorphous` | doped and hydrogenated bulks at a target density, for melt-quench |
| 4 | `04_sliding_interface` | slab / lubricant / slab with LAMMPS groups and load |
| 5 | `05_tip` | hemispherical tip, passivation, tip on a slab |

Also here: `bulk/` (the bulk cells used by the examples) and `tools/` (the
scripts, see `tools/README.md`). The orthorhombic supercell script is in
`qe/tools/ortho_cell.py` (explained in `qe/06_slab_al111`).

## Setup

The scripts need ase, numpy and, for slabs and interfaces, pymatgen: all
in the shared Python environment of `../leonardo/README.md`. The
lubricant needs packmol; on Leonardo use the group build:

```bash
export PATH=/leonardo/pub/userexternal/epedrett/software/packmol-21.2.1:$PATH
```

(add it to your `~/.bashrc`), or pass the executable with `--packmol`.
Elsewhere, build it from https://github.com/m3g/packmol (`make`, needs
gfortran).

All the examples run in seconds on a login node. Commands run from the
example folder; `-h` prints the usage of every script.

## Conventions

- **Atoms sorted along z.** Every script sorts the output along z; in
  stacked systems each part (bottom slab, molecules, top) is sorted on its
  own. LAMMPS groups are then contiguous id ranges (`group name id <> 1 72`).
- **LAMMPS data files**: `atom_style atomic`, masses included, atom types
  in the order of `pair_coeff` (`--specorder`). For a MLIP that is the
  element list of the model: write it in full (e.g. `--specorder Fe C O H`)
  even if some elements are absent, so that one `pair_coeff` line fits
  every system of the project.
- **Look at every structure** in OVITO before using it: terminations,
  overlaps, vacuum, groups (colour by id range).
- Cells: c perpendicular to the surface, vacuum along z. Sliding needs an
  orthorhombic cell with the sliding direction along x or y.
