# 01 · First SCF: bulk diamond

A self-consistent field (SCF) calculation solves the Kohn-Sham equations
for fixed atoms and gives the total energy, forces, stress and electronic
levels. Every other calculation (relaxation, NEB, single points for
MLIPs) is a series of SCF runs, so this is the place to learn how to run
pw.x and read its output. The system is diamond, 2 atoms: seconds of run
time.

The folder contains the input `C_diamond.pwi` (every line commented: read
it first), the jobscript `job.sh` (1 GPU, debug queue) and the expected
output in `reference/`.

## Run

Work in your copy of the repository (`../README.md`, "Setup"):

```bash
cd qe_runs/01_scf_diamond
qe2xyz C_diamond.pwi          # look at the structure in OVITO first
sbatch job.sh C_diamond.pwi
squeue --me                   # PD = waiting, R = running
```

It takes seconds once it starts. Output: `C_diamond.pwo`.

## Read the output

Quick status first:

```bash
python ../tools/pwo_summary.py C_diamond.pwo
```

Then open the output (`less C_diamond.pwo`, `/` to search) and find these
lines, top to bottom:

| Line in the .pwo | Expected | Meaning |
|---|---|---|
| `48 Sym. Ops., with inversion` | 48 | full cubic symmetry found: the structure is right |
| `number of k points=    29` | 29 | irreducible points of the 8×8×8 grid |
| `the Fermi energy is` | 16.05 eV | printed with smearing; in an insulator it falls somewhere in the gap and means little |
| `!    total energy` | -36.8675 Ry | comparable only with runs with the same pseudo, cutoff and k-grid |
| `smearing contrib. (-TS)` | 0.00000000 Ry | no partially occupied states: diamond has a gap |
| `convergence has been achieved in N iterations` | ~7 | SCF converged |
| `Total force` | 0.000000 | forces vanish by symmetry |
| `total   stress ... P=` | ~18 kbar | positive: PBE wants a slightly larger cell than the experimental 3.567 Å |
| `JOB DONE.` | present | the run ended; it is printed also when the SCF fails, so check the convergence line too |

Differences of ~1e-6 Ry from the reference are normal.

## Try

1. Change `A` to 3.60: what happens to the sign of `P=`?
2. Switch to `occupations = 'fixed'` (delete the `smearing` and `degauss`
   lines) and add `nbnd = 8`: the output prints `highest occupied, lowest
   unoccupied level` instead of the Fermi energy (13.35 and 17.50 eV: a gap
   of ~4.15 eV on this grid, as expected from PBE; experiment 5.5 eV). The
   total energy does not change: with a gap, smearing has no effect.
3. Add `nstep = 0` in `&CONTROL` (dry run): pw.x stops after the
   initialization. Which of the lines of the table are still there?
