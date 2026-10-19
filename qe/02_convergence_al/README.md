# 02 · Convergence tests: cutoff, k-points, smearing (fcc Al)

The plane-wave cutoff, the k-point grid and (for metals) the smearing are
numerical parameters: results must not depend on them. A convergence test
repeats the same calculation while changing one parameter, and picks the
cheapest value whose result is close enough to the converged one. Do it
for every new system and pseudopotential; here on fcc Al, 1 atom, so a
whole scan takes minutes.

The folder contains the template `Al_fcc.pwi` (only the scanned parameter
changes from run to run), `job.sh` (runs every `*.pwi` of the current
folder in one job), `kshift/` (section 5) and the reference tables and
plots in `reference/`.

## Before starting

- Converge the **quantity you need**. Energy differences (adsorption,
  barriers, surface energies) converge faster than absolute energies;
  pressure and forces are slower.
- Suggested criteria (adjust to the project): energy < 1 meV/atom;
  pressure < 0.5 kbar (needed for vc-relax and for MLIP stress labels).
- First the cutoff (at a reasonable k-grid), then the k-grid (at the chosen
  cutoff).
- Metals: k-grid and smearing go together.

## 1. Cutoff

`ecutwfc` sets the size of the plane-wave basis, `ecutrho` that of the
charge density (8-12 × `ecutwfc` for US/PAW pseudopotentials). Scan both
together first:

```bash
mkdir ecut
cd ecut
python ../../tools/make_scan.py ../Al_fcc.pwi --ecut 20 25 30 35 40 50 60 80 --dual 8
```

`make_scan.py` writes `ecut_020.pwi`, ..., each with `ecutwfc` set and
`ecutrho = 8 × ecutwfc` (`--dual 8`, for US/PAW pseudopotentials). Open
two of them and compare (`diff ecut_020.pwi ecut_080.pwi`): only two
lines differ.

```bash
sbatch ../job.sh                             # all ecut_*.pwi, one after the other
python ../../tools/conv_plot.py ecut_*.pwo   # when the job has finished
```

`conv_plot.py` prints, for each run, the energy per atom and the
difference from the last (largest) cutoff, the pressure, and saves `conv.png`.

Then the density cutoff alone, at the chosen `ecutwfc`:

```bash
mkdir ../rho
cd ../rho
python ../../tools/make_scan.py ../Al_fcc.pwi --param ecutrho 160 240 320 400 480
```

## 2. k-points

The k-grid samples the Brillouin zone; metals need dense grids. Put the
chosen cutoffs in `Al_fcc.pwi`, then:

```bash
mkdir ../k
cd ../k
python ../../tools/make_scan.py ../Al_fcc.pwi --kpts 6 8 10 12 14 16 18 20 24 28
sbatch ../job.sh
python ../../tools/conv_plot.py k_*.pwo      # when the job has finished
```

The reference is the densest grid, so the scan must go well beyond the
grid you will choose.

## 3. Smearing (optional)

In a metal the occupations are smeared over `degauss` (Ry) around the
Fermi level, so that a finite k-grid converges. The physical result is
the limit `degauss` → 0: the reference is therefore the **smallest**
value, unlike cutoff and k-grid. A larger `degauss` converges faster in k
but moves the result away from that limit; small values need denser
k-grids, so run this scan with the grid of section 2 (the reference used
16³: its smallest values may carry some k-point error).

```bash
mkdir ../degauss
cd ../degauss
python ../../tools/make_scan.py ../Al_fcc.pwi --param degauss 0.002 0.005 0.00735 0.01 0.02 0.05 0.1
sbatch ../job.sh
python ../../tools/conv_plot.py degauss_*.pwo --ref min   # when finished; reference = smallest degauss
```

## Results (reference/)

| Scan | Result |
|---|---|
| `ecutwfc` (dual 8) | within 1 meV/atom of 80 Ry from **30 Ry**; pressure within 0.25 kbar already at 25 Ry |
| `ecutrho` at 40 Ry | no effect between 160 and 480 Ry (0.005 meV/atom): for this pseudo even dual 4 is enough |
| k-grid | energy oscillates by ±20 meV/atom up to 16³; within 1 meV/atom from **20³**; pressure within 0.5 kbar from 18³ |
| `degauss` at 16³ | with respect to 0.002 Ry: −0.5 meV/atom at 0.00735 Ry (0.1 eV), −1.1 at 0.01, −5.6 at 0.02; pressure within 0.8 kbar |

Chosen for the next examples: `ecutwfc = 40`, `ecutrho = 320` (margin
over 30 Ry) and 20×20×20 k-points.

Lessons:
- The k-convergence of a metal is **not monotonic**: the energy oscillates
  around the converged value (`reference/conv_k.png`). Judge on the whole
  curve, not on the last two points.
- Whether `ecutrho` matters depends on the pseudopotential: here it does
  not; for Fe or Mn it does (SSSP uses dual 12 there).
- Al alone would allow 30 Ry, but the cutoff of a calculation is set by
  its **hardest element**: as soon as H enters (07, 08) we need 60 Ry.

## 4. Reuse the grid: k-spacing

The converged grid holds for this cell only. Its k-spacing holds for any
cell of the same material (conventional cell, supercell, slab):

```bash
python ../../tools/kgrid.py k_20.pwo --grid 20 20 20
```

prints a spacing of 0.135 Å⁻¹, used in 06-08 (`--kspacing 0.135`).
`kgrid.py` reads the `.pwo` because ASE cannot read inputs with `ibrav ≠ 0`.

## 5. Grid shift: use Γ-centered grids

The last three numbers of `K_POINTS automatic` shift the grid by half a
step (`1`) or not (`0`, the grid contains Γ). On hexagonal cells a shifted
grid silently changes the sampling. `kshift/` contains a 3-layer Al(111)
slab with one H atom away from the high-symmetry sites: the only symmetry
left is k ↔ −k (`No symmetry found` in the output). `gamma.pwi` and
`shifted.pwi` differ only in the `K_POINTS` line (`diff gamma.pwi shifted.pwi`).

1. Dry run of both (`nstep = 0` in `&CONTROL`; from `kshift/`,
   `sbatch ../job.sh` runs every `.pwi` of the folder): compare
   `number of k points` with the 7×7 = 49 points of the grid.
2. Copy both, add `nosym = .true.` in `&SYSTEM` and change `prefix`;
   dry run again.
3. Remove `nstep = 0`, delete the dry-run outputs (`rm *.pwo`: `job.sh`
   skips inputs whose `.pwo` already says `JOB DONE`), run the four inputs
   (seconds) and compare the energies.

| Input | k-points | k-points with `nosym` | E − E(`nosym`) |
|---|---|---|---|
| `gamma.pwi` (7×7×1, Γ-centered) | 25 | 25 | 0 |
| `shifted.pwi` (7×7×1, shifted) | **75** | 25 | −11.8 meV |
| same structure, 8×8×1, Γ-centered | 34 | 34 | 0 |
| same structure, 8×8×1, shifted | **96** | 32 | +10.5 meV |
| clean slab (12 symmetry operations), 7×7×1, shifted | 16 | 25 | −0.1 meV |

All values in `reference/kshift.txt`. These grids are far from converged
for Al (this cell needs 19×19×1, see section 4): only the comparison at the
same grid matters here.

Why more points than the grid itself: pw.x reduces the grid with the
symmetry operations of the **Bravais lattice** (24 for hexagonal), then
unfolds each point onto its star with the same operations and folds back
only with the operations of the **crystal**. A Γ-centered grid with
nk1 = nk2 is mapped onto itself by every operation of the hexagonal
lattice, so nothing changes. A shifted grid is not: the ±120° rotations
move its points off the grid, and the unfolding adds them. The 75 points
are the 25 of the grid plus their copies rotated by ±120°, each with a
third of the weight. With `nosym` pw.x skips this and uses only k ↔ −k.
(A Γ-centered grid with nk1 ≠ nk2 has the same problem: 8×4×1 on this cell
gives 34 points, 18 with `nosym`.)

Consequences of a shifted grid on a hexagonal cell:
- up to 3× the cost;
- the sampling depends on the symmetry of each structure and on `nosym`.
  The same structure on the same grid gives energies 11 meV apart, so a
  high-symmetry clean slab and the slab with an adsorbate, or a run with
  symmetry and a `nosym` MLIP single point, are not sampled the same way
  and the error does not cancel.

With Γ-centered grids symmetry changes only the cost: same energy with and
without `nosym`. VASP refuses shifted grids on hexagonal lattices (see
[KPOINTS](https://vasp.at/wiki/KPOINTS), "Symmetry reduction of the
mesh"); pw.x accepts them. Use `0 0 0`: all scripts in `tools/` do.

## Exercise: a magnetic metal

Repeat the cutoff scan with `04_lattice_gallery/bcc_Fe.pwi` (`nspin = 2`).
Compare the cutoff you need with Al, check `ecutrho` separately, and check
the magnetization in each run (`python ../../tools/pwo_summary.py *.pwo`,
column `M_tot`, ~2.2-2.3 μB/atom).

## Common problems

- `SCF-FAIL` at very low cutoff or coarse k: drop that point, or lower
  `mixing_beta` to 0.3.
- A run killed by the 30 min limit: resubmit `sbatch ../job.sh`; finished
  runs are skipped.
- All runs give the same energy: the scanned variable is not in the
  template, or appears twice (`make_scan.py` stops in that case).
