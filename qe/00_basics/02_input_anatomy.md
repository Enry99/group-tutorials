# Anatomy of a pw.x input

A pw.x input is a text file of namelists (variables, `name = value`) and
cards (atoms, positions, k-points). This page gives the layout, the units,
a checklist before submitting, the variables you change most often and a
few pw.x behaviours worth knowing. Reference for every variable:
[pw.x guide](https://www.quantum-espresso.org/Doc/INPUT_PW.html).

## Structure

```
&CONTROL ... /        what to compute, where files go
&SYSTEM ... /         cell, atoms, cutoffs, occupations, spin, vdW
&ELECTRONS ... /      SCF algorithm and threshold
&IONS ... /           needed for relax, md, vc-relax
&CELL ... /           needed for vc-relax, vc-md
ATOMIC_SPECIES        label, mass, pseudopotential file
CELL_PARAMETERS       only with ibrav = 0
K_POINTS              k-point grid
HUBBARD               only for DFT+U
ATOMIC_POSITIONS      positions (+ optional 0/1 flags to fix coordinates)
```

Namelists (`&NAME ... /`) come in this order. Cards (the uppercase blocks)
come after them, in any order; in these examples `ATOMIC_POSITIONS` is
always last, so the positions are at the end of the file.

**Units.** Energies in Ry, lengths in bohr, unless stated otherwise.
1 Ry = 13.6057 eV; 1 bohr = 0.529177 Å; force 1 Ry/bohr = 25.711 eV/Å;
stress in the output in kbar (1 GPa = 10 kbar). `A`, `C` and
`ATOMIC_POSITIONS angstrom` are in Å; `celldm(1)` is in bohr.

<!-- ## Comments (tested with pw.x 7.3.1)

| Where | Works? |
|---|---|
| `!` after a value in a namelist, or on its own line inside a namelist | yes |
| `!` at the end of a data line of a card (`Fe 0.0 0.0 0.0 ! fixed`) | yes |
| comment line **between** two cards | yes |
| comment line **inside** a card, between its data lines | **no**: `cannot read atomic specie from: ! ...` or `Error while parsing atomic position card.` |
| blank lines between cards or inside a card | yes | -->

## Checklist before submitting

1. **Look at the structure** every time you build or edit one: most
   "strange" errors are wrong structures. `qe2xyz file.pwi` writes an
   `.extxyz` for OVITO (`01_running_qe.md`, "Looking at structures"); it also
   warns about atoms closer than 0.7 Å.
2. **Dry run** for new or expensive inputs: with `nstep = 0` in `&CONTROL`
   pw.x stops right after the initialization (seconds). It reads the
   pseudopotentials and prints cell, symmetry operations and number of
   k-points. Submit it with the 1-GPU debug jobscript (the GPU build may
   not start on login nodes, which have no GPU), check, then set the real
   `nstep` and submit.
3. `nat`, `ntyp` match the cards; one pseudo per species; same pseudos
   and cutoffs as the runs you will compare with.
4. Magnetic? `nspin = 2` and sensible `starting_magnetization`.
5. `max_seconds` ~5% below the walltime; `nstep` set for relax/md;
   `disk_io` chosen.

## Variables you will touch in every run

**&CONTROL**

| Variable | Meaning |
|---|---|
| `calculation` | `'scf'`, `'relax'` (atoms), `'vc-relax'` (atoms + cell), `'md'`, `'nscf'`, `'bands'` |
| `prefix`, `outdir` | name and folder of the files pw.x writes (wavefunctions, charge, restart). Two runs in the same folder need different `prefix` or `outdir` |
| `pseudo_dir` | folder with the `.UPF` files (the group folder is already set in the inputs) |
| `tprnfor`, `tstress` | print forces / stress (automatic for relax and vc-relax) |
| `max_seconds` | stop cleanly before the walltime: always set it, ~5% below the walltime of the jobscript |
| `disk_io` | what pw.x writes to `outdir` (table below) |
| `forc_conv_thr`, `etot_conv_thr` | relax thresholds (Ry/bohr, Ry) |
| `nstep` | maximum ionic steps. **Default 50 for relax, vc-relax and md** (1 for scf): often not enough. Always set it, e.g. `nstep = 200`. |

`disk_io` (QE 7.3; default `'low'` for scf, `'medium'` otherwise):

| Value | Writes | Use |
|---|---|---|
| `'none'` | nothing | scf/single points whose wavefunctions and density you do not need (MLIP labels, convergence scans, E(V)). No restart |
| `'nowf'` | xml + charge density, no wavefunctions | when you need only the density (e.g. to restart the SCF of another run from it). No restart of an interrupted run |
| `'low'` | charge at each SCF step, wavefunctions only at the end | relaxations in many folders: less disk than `'medium'`, restart still possible |
| `'medium'` | as `'low'`, wavefunctions on disk during the run if a process has several k-points | default for relax/md |
| `'high'` | everything at every step | almost never needed |

Wavefunction files of a slab with many k-points take GBs: delete `outdir`
(`WORK/`) when a run is finished and you do not need to restart it.

**&SYSTEM**

| Variable | Meaning |
|---|---|
| `ibrav`, `A`/`celldm(1)`, `C`/`celldm(3)`, ... | lattice (see `04_lattice_gallery`) |
| `nat`, `ntyp` | number of atoms and of species |
| `ecutwfc`, `ecutrho` | cutoffs (Ry); `ecutrho` = 8-12 × `ecutwfc` for US/PAW, 4 × for NC |
| `occupations` | `'fixed'` (insulators, isolated atoms and molecules; with `nspin = 2` it needs `tot_magnetization`), `'smearing'` (metals, magnetic systems, anything unknown) |
| `smearing`, `degauss` | smearing function and width (Ry) |
| `nspin`, `starting_magnetization(i)` | spin-polarized run; initial moment of species i |
| `vdw_corr` | dispersion correction, e.g. `'dft-d3'` |
| `nbnd` | number of bands (sometimes useful to increase for Fe or complex systems, may improve convergence) |
| `nosym` | do not use symmetry (MLIP single points, relax for adsorption, MD) |

**&ELECTRONS**

| Variable | Meaning |
|---|---|
| `conv_thr` | SCF threshold on the **total** energy (Ry): it is extensive, scale it with the number of atoms |
| `mixing_beta` | 0.7 default; 0.1-0.3 for magnetic metals, oxides, slabs |
| `mixing_mode` | `'plain'` (bulk), `'local-TF'` (slabs, interfaces, molecules on surfaces) |
| `mixing_ndim` | mixing history, default 8; 16-20 suggested |
| `electron_maxstep` | maximum SCF iterations (default 100, increase it to at least 200 to be safe) |

**K_POINTS**

```
K_POINTS automatic
8 8 8  0 0 0          ! nk1 nk2 nk3, then the shifts: 0 = Γ-centered, 1 = shifted
```

- Use Γ-centered grids (`0 0 0`), with nk1 = nk2 on hexagonal cells
  (fcc(111) and hcp(0001) slabs, graphite, MoS₂, MXenes, ...). Otherwise
  pw.x uses up to 3× more k-points, and the sampling depends on the
  symmetry: see `02_convergence_al`, section 5. The scripts in `tools/`
  write Γ-centered grids.
- Slabs: `nk3 = 1`, no shift along z.
- `K_POINTS gamma`: Γ only, with a faster algorithm (real wavefunctions);
  for large cells and molecules. Not the same as `1 1 1 0 0 0` (same
  point, slower algorithm).
- Grids from VASP (e.g. Materials Project): `Gamma N` = QE `N 0`;
  `Monkhorst-Pack N` = QE `N 0` for odd N, `N 1` for even N.


## Behaviours worth knowing (tested with pw.x 7.3.1)

- A variable in the wrong namelist, or a typo, stops pw.x:
  `bad line in namelist &control: "  ecutwfc = 30" (error could be in the previous line)`.
- A variable written twice is accepted silently; the **last** value wins.
- `calculation = 'relax'` without `&IONS` (or vc-relax without `&CELL`):
  `bad line in namelist &ions:` followed by a line of a card. The error
  points at the card because pw.x went on looking for the namelist.
- `nspin = 2` without `starting_magnetization`: with smearing,
  `some starting_magnetization MUST be set`; with fixed occupations,
  `fixed occupations and lsda need tot_magnetization`.
- A metal with `occupations = 'fixed'` runs **without error** and gives a
  wrong result. Check the output: metals print `the Fermi energy is`,
  fixed occupations print `highest occupied level`.
- `starting_magnetization`: values in (-1, 1) are fractions of the
  valence charge. If any value has |value| ≥ 1, pw.x 7.3 divides all of
  them by the valence charge (moments in μB). Use one convention per input.
- `ecutrho < 4*ecutwfc` only prints `ecutrho < 4*ecutwfc, are you sure?` and
  continues.
- `JOB DONE.` is printed also when the SCF did **not** converge
  (`convergence NOT achieved after N iterations: stopping` a few lines
  above). Always check convergence, e.g. with `tools/pwo_summary.py`.
- Old DFT+U variables (`lda_plus_u`, `Hubbard_U`) stop pw.x since 7.1:
  `DFT+Hubbard input syntax has changed since v7.1`; use the `HUBBARD` card.
- With `tefield`/`dipfield` (dipole correction) pw.x prints `Presently no
  symmetry can be used with electric field`: no symmetry, more k-points.
- Errors are also written to a file `CRASH` in the run folder.
