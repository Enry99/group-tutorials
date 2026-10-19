# 04 · Gallery of common crystal structures

One ready-to-run SCF input per structure type. Copy the one closest to your
material, change element, lattice parameters and internal coordinates, then
converge cutoff and k-points (02) and the lattice (03).

All files were checked: pw.x 7.3.1 reads them and finds the expected
symmetry, spglib gives the space group below, and the nearest-neighbour
distances match the experimental structures. fcc Al is in 02.

| File | Prototype | Space group | ibrav | nat | Check (NN distance) |
|---|---|---|---|---|---|
| `fcc_Ni.pwi` | A1 fcc | Fm-3m (225) | 2 | 1 | Ni-Ni 2.49 Å |
| `bcc_Fe.pwi` | A2 bcc | Im-3m (229) | 3 | 1 | Fe-Fe 2.48 Å |
| `bcc_Fe_conventional.pwi` | A2, cubic cell | Im-3m (229) | 0 | 2 | same energy/atom as `bcc_Fe` |
| `hcp_Ti.pwi` | A3 hcp | P6₃/mmc (194) | 4 | 2 | Ti-Ti 2.90 Å |
| `diamond_Si.pwi` | A4 diamond | Fd-3m (227) | 2 | 2 | Si-Si 2.35 Å |
| `zincblende_SiC.pwi` | B3 | F-43m (216) | 2 | 2 | Si-C 1.89 Å |
| `rocksalt_TiC.pwi` | B1 | Fm-3m (225) | 2 | 2 | Ti-C 2.17 Å |
| `rutile_TiO2.pwi` | C4 | P4₂/mnm (136) | 6 | 6 | Ti-O 1.95 Å |
| `corundum_Al2O3.pwi` | D5₁, hexagonal cell | R-3c (167) | 4 | 30 | Al-O 1.85 Å |
| `wurtzite_ZnO.pwi` | B4 | P6₃mc (186) | 4 | 4 | Zn-O 1.97 Å |
| `layered_MoS2_2H.pwi` | C7, 2 layers | P6₃/mmc (194) | 4 | 6 | Mo-S 2.42 Å |
| `layered_graphite_AB.pwi` | Bernal graphite | P6₃/mmc (194) | 4 | 4 | C-C 1.42 Å |

## How QE describes a crystal

**Bravais lattice.** `ibrav` selects the lattice and its orientation; the
lattice parameters are given either in Å (`A`, `B`, `C`, `cosAB`, `cosAC`,
`cosBC`) or in bohr and ratios (`celldm(1)` = a in bohr, `celldm(3)` = c/a,
...). Never both: pw.x stops with `do not specify both celldm and a,b,c!`.
The primitive vectors for each `ibrav` are listed in
[INPUT_PW](https://www.quantum-espresso.org/Doc/INPUT_PW.html#ibrav): read
them before writing crystal coordinates, they are not always the ones of
Materials Project or of a textbook.

| ibrav | Lattice | Parameters |
|---|---|---|
| 1 / 2 / 3 | simple / face-centered / body-centered cubic | `A` (conventional cube edge) |
| 4 | hexagonal (also used for trigonal crystals in the hexagonal setting) | `A`, `C` |
| 5 | trigonal R, rhombohedral primitive cell | `A`, `cosAB` |
| 6 / 7 | simple / body-centered tetragonal | `A`, `C` |
| 8 | simple orthorhombic | `A`, `B`, `C` |
| 0 | any cell, given in `CELL_PARAMETERS` | vectors in `angstrom`, `bohr` or `alat` |

**Positions.** `ATOMIC_POSITIONS crystal` (fractions of the cell vectors) is
the safest choice: it does not change when you scale the cell. `angstrom`
and `bohr` are Cartesian; `alat` is Cartesian in units of the lattice
parameter (the course used it).

**ibrav ≠ 0 or ibrav = 0?**
- ibrav ≠ 0: readable, and vc-relax with `cell_dofree = 'ibrav'` keeps the
  lattice type.
- ibrav = 0: what ASE writes and what you need for supercells, slabs,
  interfaces and anything built by a script. pw.x finds the symmetry anyway.
- **ASE cannot read a `.pwi` with ibrav ≠ 0** (`ValueError: ASE does not
  support ibrav != 0`). Use `tools/qe2xyz.py`, or read the `.pwo` (it
  contains the cell).

## Notes per structure

- **Magnetic metals** (Fe, Ni, Co): `nspin = 2` and a `starting_magnetization`,
  otherwise the result is non-magnetic and wrong.
- **Primitive vs conventional cell**: `bcc_Fe_conventional.pwi` has 2 atoms
  instead of 1, and the same energy per atom when the k-point *spacing* is
  similar (checked: 0.9 meV/atom with 16³ vs 12³). Use the conventional cell
  only when you need it (cubic supercells, orthogonal slabs).
- **Corundum** in the hexagonal cell is three primitive cells: pw.x prints
  `This is a supercell, fractional translations are disabled` and uses 6
  symmetry operations instead of 12. Fine for building (0001) slabs; for
  bulk properties the 10-atom rhombohedral cell (ibrav = 5) is cheaper.
  α-Fe₂O₃ has the same structure plus antiferromagnetic order: 05.
- **Layered** MoS₂ and graphite need `vdw_corr = 'dft-d3'`: plain PBE
  barely binds the layers.
- **Internal coordinates** (u in rutile, wurtzite, corundum, z of S in MoS₂)
  are experimental: relax them (`relax` or `vc-relax`) before using the
  structure.
- **k-grids** here correspond to spacings of ~0.18-0.20 Å⁻¹ for metals and
  ~0.25-0.30 Å⁻¹ for insulators (`tools/kgrid.py`): starting points, not
  converged values. Al needed 0.135 Å⁻¹ for 1 meV/atom (02).
- **Pseudopotentials**: group set where it exists, otherwise pslibrary 1.0.0;
  see `00_basics/03_pseudopotentials.md`.

## Check any structure before running

```bash
qe2xyz rutile_TiO2.pwi                              # -> .extxyz, open in OVITO (any ibrav)
python ../tools/struct_to_qe.py my.cif -o my.pwi    # space group, primitive or not
```

`qe2xyz` (alias of `tools/qe2xyz.py`, see `00_basics/01_running_qe.md`)
builds the cell with QE's own formulas: checked against pw.x for every
ibrav from 1 to 14 (including the negative variants and 91) and for all
files of this folder.
