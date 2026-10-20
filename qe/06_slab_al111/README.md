# 06 · Surface slab: Al(111), relaxation and surface energy

A surface is modelled by a slab: a few atomic layers of the crystal,
periodic in the plane, separated from their periodic images along z by
vacuum. The slab is cut from the relaxed bulk, relaxed, and its energy
compared with the bulk gives the surface energy γ, the cost of creating the
surface. γ and the relaxations depend on the number of layers: the last
step converges them with the thickness. The same slab is the starting point
for adsorption (07) and sliding interfaces.

The folder contains the 6-layer slab `Al111_6L.pwi` written by
`tools/build_slab.py` (a0 = 4.012 Å from 03), `job.sh` (1 GPU, debug
queue) and the surface energies for 3 to 6 layers in `reference/`. All
atoms relax (no fixed layers), and the atoms are sorted along z.

## 1. Build

`tools/build_slab.py` cuts the slab with ASE and writes a pw.x relax input
with the settings of the bulk:

```bash
python ../tools/build_slab.py Al fcc111 --a0 4.012 --layers 6 \
       --template ../03_lattice_parameter/01_eos/Al_fcc.pwi --kspacing 0.135
```

What the arguments do:

- `Al fcc111`: element and surface (ASE builders: also fcc100, fcc110,
  bcc100, bcc110, bcc111, hcp0001).
- `--a0 4.012`: the lattice parameter from 03, computed with the same
  settings.
- `--template`: the bulk input. Cutoffs, smearing, D3 and pseudopotential
  are copied from it, so slab and bulk are consistent.
- `--kspacing 0.135`: the bulk k-spacing from 02. The script turns it into
  an in-plane grid for this cell (19×19×1) with one k-point along z.

It also sets: 15 Å of vacuum, `calculation = 'relax'`, `nstep = 200`,
`forc_conv_thr = 1.0d-3` Ry/bohr (≈ 0.026 eV/Å), `mixing_mode = 'local-TF'`
(better for metal + vacuum), `conv_thr` scaled with the number of atoms,
and writes `slab_info.json` with all the parameters. The first line of the
`.pwi` records the command that wrote it. Open the `.pwi` and check each of
these lines.

**Look at the slab before running:**

```bash
qe2xyz Al111_6L.pwi        # -> Al111_6L.extxyz, open it in OVITO
```

## 2. Relax

```bash
sbatch job.sh Al111_6L.pwi
python ../tools/pwo_summary.py Al111_6L.pwo     # when the job has finished: DONE?
```

## 3. Surface energy

γ = (E_slab − N e_bulk) / (2A): the energy cost per area of the two
surfaces created. e_bulk must come from the same settings. The final scf of
the vc-relax in 03 is at the relaxed a0 and has the same cutoff, smearing
and D3:

```bash
python ../tools/surface_energy.py Al111_6L.pwo \
       ../03_lattice_parameter/02_vc_relax/Al_fcc_vcrelax.pwo
```

Both surfaces of an fcc(111) slab are equivalent and both relax, so the
formula is exact here.

## 4. Converge with the thickness

Build and relax slabs with 3 to 6 layers (one folder each), then compute γ
for each. Reference results (`reference/gamma_vs_layers.txt`):

| Layers | γ (J/m²) | Δd₁₂ |
|---|---|---|
| 3 | 0.973 | +0.2% |
| 4 | 1.041 | +2.0% |
| 5 | 1.024 | +2.8% |
| 6 | 1.024 | +1.6% |

- γ is converged within 0.01 J/m² from 5 layers. Thin slabs oscillate
  (3 layers is 5% low, 4 layers 2% high): judge on at least three
  thicknesses, not two.
- Δd₁₂, the change of the first interlayer distance, is positive: Al(111)
  is a known exception, most metal surfaces contract. It converges more
  slowly than γ: relaxations need thicker slabs than surface energies.
- Check e_bulk: E(N) − E(N−1) must tend to e_bulk. Here E(6) − E(5) =
  −537.6099 eV and e_bulk = −537.6105 eV.
- D3 matters: without D3 (a0 = 4.039 Å, 3 layers) γ = 0.78 J/m². Compare
  with literature only at the same functional and D3 setting.

## Orthorhombic cells (sliding MD, interfaces)

Hexagonal and bcc(110) surface cells are oblique (60/120°, 54.7°). Sliding
simulations and many interface builders need an orthorhombic cell with the
sliding direction along x or y. `tools/ortho_cell.py` rotates a1 onto x,
builds the smallest orthorhombic supercell (2× the atoms here), repeats it
in-plane and sorts the atoms along z:

```bash
python ../tools/ortho_cell.py Al111_6L.pwo   # relaxed slab -> Al111_6L_ortho.extxyz
python ../tools/ortho_cell.py Al111_6L.pwo --square 20      # both sides close to 20 A
python ../tools/ortho_cell.py Al111_6L.pwo --repeat 4 2 -f lammps-data -o slab.data
```

It checks that the cell is orthogonal, that the number of atoms is the
expected one and that no atom was duplicated. For a QE input of the new
cell, pass the `.extxyz` to `tools/struct_to_qe.py`.

The script prints the supercell matrix P. ASE and pymatgen convention: the
rows of P are the new lattice vectors in units of the old ones,
A_new = P·A, with the lattice vectors as rows of A (so P = A_new·A⁻¹).
VESTA uses the transpose (lattice vectors as columns, A_new = A·P).
Example, √3×√3R30° in the ASE convention: `[[2,1,0],[-1,1,0],[0,0,1]]` for a
120° cell (Materials Project, pymatgen), `[[1,1,0],[-1,2,0],[0,0,1]]` for
the 60° cells of the ASE slabs used here. The matrix depends on the setting
of the cell: check the result (cell angles, number of atoms) every time.

## Other surfaces

`build_slab.py` covers elemental fcc, bcc and hcp surfaces (the ASE
builders). Compounds, other Miller indices and the choice among several
terminations need a general slab generator (pymatgen `SlabGenerator`):
see `structures/01_slabs`.

## Common problems


- γ does not converge with the thickness but drifts linearly with N: e_bulk
  from a run with a different cutoff, smearing or k-spacing. The error
  grows with N instead of cancelling.
- SCF slow or not converging: `mixing_beta` 0.1-0.2; magnetic slabs need
  more patience (and a k-grid that is not too coarse).
- The relax stops with `Maximum CPU time exceeded` (thick slabs in the
  30-min debug queue): continue it with `restart_mode = 'restart'`
  (`troubleshooting.md`, "Time limits and restart"). The BFGS history is
  kept: the step count continues from where it stopped.
- Many BFGS steps with tiny energy changes: normal for metal slabs, whose
  energy surface is flat (12-16 steps here).
