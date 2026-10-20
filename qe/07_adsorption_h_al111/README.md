# 07 · Adsorption energy: H on Al(111)

The adsorption energy says how strongly a species X binds to a surface:

    E_ads = E(slab+X) − E(slab) − E(X)

X is the adsorbate as an isolated species: an atom or a molecule alone in
a large box. E_ads < 0 means that X is more stable on the surface than in
the gas phase. The three energies come from three calculations that must
use identical settings. Here X is one H atom on Al(111), on two sites (the
fcc and hcp hollows); the exercise at the end adsorbs a whole H₂ molecule,
dissociated on the surface.

The work is split into three folders:

- `01_gas_phase`: the isolated H atom (`H_atom.pwi`), and the H₂ molecule
  for the exercise (`H2_box.pwi`);
- `02_clean_slab`: the Al(111) slab without H, relaxed;
- `03_slab_h`: H on the relaxed slab, relaxed on each site.

01 and 02 are independent: submit them together. 03 starts from the
relaxed slab of 02. The commands of each section start from this folder.

## The cutoff changes here

Al alone needs 30 Ry (02), H needs more (its pslibrary PAW suggests
46 Ry; SSSP uses 60 Ry). The hardest element sets the cutoff of every
calculation that enters E_ads, so all three steps use 60/480 Ry, including
the clean slab: it cannot be taken from 06, which used 40 Ry.

## 1. Gas phase: the isolated H atom

One atom: a single scf in a 12 Å box, Γ point only. The H atom has one
unpaired electron, so the run is spin-polarized with the moment fixed to
1 μB (`nspin = 2`, `tot_magnetization = 1`). Without spin polarization its
energy is 1.1 eV too high, and every E_ads computed from it is off by
the same amount. Two more differences from the slab runs: fixed
occupations (an isolated atom has a gap, smearing is not needed) and a
looser `conv_thr` (1e-6 Ry is negligible for one atom, and much tighter
thresholds can stall on the numerical noise of the nearly empty box).

```bash
cd 01_gas_phase
sbatch job.sh H_atom.pwi
grep "total magnetization" H_atom.pwo | tail -1    # when finished: 1.00 Bohr mag/cell
```

Reference output in `reference/` (E = −1.003964 Ry).

For a molecule the procedure is the same, with a relaxation in the box
(`H2_box.pwi`; reference output in `reference/`, H-H = 0.750 Å).

## 2. The clean slab

```bash
cd 02_clean_slab
python ../../tools/build_slab.py Al fcc111 --a0 4.012 --layers 4 --repeat 2 2 \
       --ecut 60 480 --template ../../03_lattice_parameter/01_eos/Al_fcc.pwi --kspacing 0.135
```

Compared with 06: `--repeat 2 2` (a 2×2 cell, so H atoms of neighbouring
cells are 5.7 Å apart, i.e. 1/4 monolayer) and `--ecut 60 480`. The k-grid
becomes 10×10×1 for the same k-spacing.

```bash
sbatch job.sh Al111_4L_2x2.pwi
```

## 3. H on two sites

When step 2 is finished (`python ../tools/pwo_summary.py 02_clean_slab/Al111_4L_2x2.pwo`
says `DONE`):

```bash
cd 03_slab_h
python ../../tools/add_adsorbate.py ../02_clean_slab/Al111_4L_2x2.pwo H fcc --height 1.0
python ../../tools/add_adsorbate.py ../02_clean_slab/Al111_4L_2x2.pwo H hcp --height 1.0
```

`add_adsorbate.py` reads the relaxed slab, puts H 1.0 Å above the top
layer on the chosen site (ASE names for fcc(111): `ontop`, `bridge`, `fcc`,
`hcp`), copies all settings of the clean-slab input and sorts the atoms
along z (H last). Look at both structures (`qe2xyz *.pwi`), then:

```bash
sbatch job.sh Al111_4L_2x2_H_fcc.pwi
sbatch job.sh Al111_4L_2x2_H_hcp.pwi
```

The two runs share the folder: they have different `prefix`, so their
files in `WORK/` do not mix.

## 4. Adsorption energy

```bash
python ../tools/ads_energy.py 03_slab_h/Al111_4L_2x2_H_fcc.pwo \
       02_clean_slab/Al111_4L_2x2.pwo 01_gas_phase/H_atom.pwo
```

The script prints the three energies with their cutoffs, warns if the
cutoffs differ or if the atoms do not add up (N(slab+X) − N(slab) = N(X)),
and gives E_ads. Repeat for `hcp`.

### Local check (reduced settings)

The chain was run locally on 2 CPU cores at reduced settings: 3 layers,
40/320 Ry and 4×4×1 k-points for all three runs. These numbers check the
procedure; they are **not** the result at the settings of this README:

| Site | E_ads | H–Al (3 neighbours) | H above the top layer |
|---|---|---|---|
| fcc | −1.91 eV | 1.91 Å | 0.88 Å |
| hcp | −1.67 eV | 1.95 Å | 1.04 Å |

Site differences need converged k-points: at the relaxed geometries,
E(hcp) − E(fcc) = 0.23 eV with 4×4×1, 0.18 eV with 8×8×1 and 0.19 eV with
10×10×1 (the k-spacing of 02).

## Consistency rules (the usual source of wrong numbers)

The three runs must share pseudopotentials, cutoffs and `vdw_corr`; slab
and slab+H also the cell, the k-grid, the smearing and the vacuum.
`add_adsorbate.py` copies them from the clean-slab input; `ads_energy.py`
checks the cutoffs and the atom count.

## Things to know

- **Reference**: E_ads depends on it. With respect to the H atom it is
  −1.9 eV (fcc, local check); with respect to ½ H₂ it is less negative by
  half the H₂ binding energy (2.2 eV at those settings) and becomes
  positive, +0.27 eV: dissociating H₂ on Al(111) costs energy.
  Always say which reference you use.
- **Dipole correction**: with H on one side only, the slab is asymmetric
  and has a small dipole. For H on a metal it is negligible and is not used
  here. For polar adsorbates or oxide surfaces add `tefield = .true.`,
  `dipfield = .true.` (`&CONTROL`) and `edir`, `emaxpos`, `eopreg`, `eamp`
  (`&SYSTEM`; see INPUT_PW), identical in slab and slab+adsorbate. Note
  that pw.x then uses no symmetry (`Presently no symmetry can be used with
  electric field`): more k-points.
- **Final site**: a relaxation started on a bridge may roll into a hollow.
  Always check the final position, not the starting label.
- **High-symmetry sites** keep their symmetry during a relaxation (the
  forces parallel to the surface are zero by symmetry): a relaxed site may
  be a saddle point or a maximum. Displace the atom by 0.1 Å and relax again
  to check. Here it happens with hcp: in the local check a NEB fcc → hcp
  shows H sliding from hcp to fcc without a barrier (08).

## Exercise: dissociative adsorption of H₂

Now X is the whole H₂ molecule, adsorbed as two H atoms:

    E_ads = E(slab+2H) − E(slab) − E(H₂)

1. Relax H₂ in its box: `sbatch job.sh H2_box.pwi` in `01_gas_phase`.
2. Add a second H to the relaxed slab+H. `add_adsorbate.py` needs
   `slab_info.json` next to the slab it reads:
   ```bash
   cd 03_slab_h
   cp ../02_clean_slab/slab_info.json .
   python ../../tools/add_adsorbate.py Al111_4L_2x2_H_fcc.pwo H fcc --offset 1 1 --height 1.0
   sbatch job.sh Al111_4L_2x2_H_fcc_H_fcc_11.pwi
   ```
   `--offset 1 1` takes the fcc site one surface cell further along a1 and
   a2: the two H are 2.84 Å apart (1/2 monolayer).
3. From `07_adsorption_h_al111`:
   ```bash
   python ../tools/ads_energy.py 03_slab_h/Al111_4L_2x2_H_fcc_H_fcc_11.pwo \
          02_clean_slab/Al111_4L_2x2.pwo 01_gas_phase/H2_box.pwo
   ```

Compare with 2 × E_ads(H) + D(H₂), where D(H₂) = 2 E(H) − E(H₂) is the
binding energy of H₂: the difference is the interaction between the two
adsorbed H atoms.

## Automating the chain (later)

When you know each step, Slurm can start step 3 after step 2 with
`sbatch --dependency=afterok:<jobid of step 2> ...`. The input of step 3
must exist by then, so it needs a small serial job that runs
`add_adsorbate.py` in between. Not needed for learning: do it by hand.
