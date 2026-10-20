# 05 · Antiferromagnetic oxide: α-Fe₂O₃ (hematite), PBE and PBE+U

Hematite is an antiferromagnetic insulator: the Fe moments alternate up
and down from layer to layer. A magnetic calculation must be given the
intended order, otherwise the SCF can end in another state (ferromagnetic,
or non-magnetic). PBE also underestimates the localization of the Fe 3d
electrons and the gap, which DFT+U partly corrects. Here: how to set up a
collinear antiferromagnet, how to check that the run ended in the intended
state, and how to add U.

Two runs of the same structure, in two folders: `01_pbe/Fe2O3_afm.pwi`
(PBE+D3, AFM order set through two Fe species) and
`02_pbe_u/Fe2O3_afm_U.pwi` (the same plus U on the Fe 3d states, `HUBBARD`
card).

The group's production choice (with or without U) depends on the project:
Fe/FeOx interfaces and MLIP datasets that also contain metallic Fe usually
avoid U, so that metal and oxide are described by the same functional.

## How the AFM order is written

pw.x has one `starting_magnetization` per **species**, not per atom. To have
spin-up and spin-down Fe, declare two species with the same pseudopotential:

```
ATOMIC_SPECIES
Fe1  55.845  Fe.pbe-sp-van.UPF
Fe2  55.845  Fe.pbe-sp-van.UPF
...
starting_magnetization(1) =  0.5    ! Fe1 up
starting_magnetization(2) = -0.5    ! Fe2 down
```

Which Fe is up: in hematite the Fe atoms form buckled layers perpendicular
to c (6 layers in the hexagonal cell, at z ≈ 0, 1/6, ..., 5/6). Spins are
parallel within a layer and alternate from layer to layer. In the input
each Fe line ends with its layer index and spin. Checked on the structure:

| Fe-Fe distance | Connection | Spins |
|---|---|---|
| 2.89 Å | face-sharing pair along c | antiparallel |
| 3.0 Å | edge-sharing, same layer | parallel |
| 3.4 Å, 3.7 Å | corner-sharing | antiparallel |

This is the experimental ground state (in the 4-Fe rhombohedral cell it is
the "+ − − +" sequence along the c axis).

**Pitfall**: a single Fe species with a positive
`starting_magnetization` starts **ferromagnetic**, and the SCF can end in a
ferromagnetic or ferrimagnetic state, higher in energy and with the wrong
electronic structure. pw.x gives no warning: the only check is the
magnetization in the output.

## Run

```bash
cd 01_pbe
qe2xyz Fe2O3_afm.pwi          # in OVITO, colour by qe_type: Fe1 and Fe2 layers alternate
sbatch job.sh Fe2O3_afm.pwi
cd ../02_pbe_u
sbatch job.sh Fe2O3_afm_U.pwi
```

`job.sh` uses one full node (4 GPUs, `-nk 4`): 30 atoms, 6×6×2 k-points.

## Check the magnetic state (always)

```bash
grep -E "total magnetization|absolute magnetization" Fe2O3_afm.pwo | tail -2
grep -A31 "Magnetic moment per site" Fe2O3_afm.pwo | tail -31
```

- `total magnetization` ≈ 0 μB/cell: the up and down moments cancel.
- `absolute magnetization` ≈ 12 × |m(Fe)| (plus a bit from O).
- In `Magnetic moment per site` the 12 Fe have the same |m| with the sign
  of their species (Fe1 +, Fe2 −); O ≈ 0.
- PBE+U: larger |m(Fe)| and a larger gap than PBE (local test values of
  the moments in "Test runs" below). Experiment: ~4.6-4.9 μB per Fe, gap
  ~2 eV; PBE underestimates both.

If the total magnetization is not ~0, the run found another magnetic state:
check the species/spin assignment, start again with larger initial moments
(±0.6-0.8), lower `mixing_beta`.

## About starting_magnetization (QE 7.3)

- Values strictly between -1 and 1: fraction of the valence electrons of
  that species (the classic meaning; used here).
- If **any** value has |value| ≥ 1, pw.x 7.3 divides **all** of them by the
  valence charge, i.e. treats them as moments in μB. Do not mix the two
  conventions in one input. ASE ≥ 3.23 writes `atoms.set_initial_magnetic_moments`
  values in μB this way (useful for slabs: different moments → separate
  `Fe`, `Fe1` species automatically).

## DFT+U notes

- Syntax since QE 7.1: a `HUBBARD {projector}` card. The old `lda_plus_u`
  / `Hubbard_U(i)` variables stop pw.x with
  `DFT+Hubbard input syntax has changed since v7.1`.
- `ortho-atomic` projectors (recommended) support forces and stress, so
  relax and vc-relax work. `atomic` is the old default.
- U must be given for **each species** that represents Fe (Fe1 and Fe2).
- U depends on the projectors and pseudopotential: values from VASP or
  Materials Project (U = 5.3 eV for Fe there) are not transferable. Either
  take a value obtained with the same projectors (4.2 eV in the input, from
  the paper cited there) or compute it (hp.x, linear response).
- Energies with different U, or with and without U, cannot be compared or
  mixed in the same dataset.

## Test runs

Local check at reduced settings, not these inputs: QE 7.3.1 (CPU), a
softer pslibrary Fe pseudo, 40/320 Ry and 2x2x1 k-points (qualitative
only):

| Run | SCF iterations | total / absolute magnetization (μB/cell) | m(Fe) in the output sphere |
|---|---|---|---|
| PBE, AFM | 62 | 0.00 / 45.9 | ±3.41 μB |
| PBE+U (U = 4.0 eV; the input now uses 4.2), AFM | 44 | 0.00 / 51.3 | ±3.80 μB |

The signs follow the Fe1/Fe2 species of the input; O carries no moment.
The per-site values are integrated in small spheres (`R=0.159` alat in the
output), so they are lower than the full moment: compare them only between
runs with the same setup. With the mixing of the inputs (`mixing_beta = 0.1`,
`mixing_ndim = 16`) the SCF needs 40-60 iterations: normal for this system.
A first attempt with a hard Fe pseudo at a cutoff far below its minimum
(30 Ry) stalled around 1e-2 Ry: an SCF that does not converge is often a
setup problem (cutoff, structure, magnetic start) before being a mixing
problem.
