# Quantum ESPRESSO quick guide

Small, runnable examples of common DFT calculations: commented inputs,
Leonardo jobscripts, what to check in the output, and the errors you will
meet. Not a QE course: for the full documentation see
[pw.x](https://www.quantum-espresso.org/Doc/INPUT_PW.html) and the
[User guide](https://www.quantum-espresso.org/Doc/user_guide/user_guide.html).

(Tested for QE version **7.3.1**).

## Contents

Each folder is one type of calculation, with a `README.md` that explains
it, the inputs, a `job.sh` and, where useful, reference results.

| # | Folder | You learn |
|---|---|---|
| 0 | `00_basics/` | running QE on Leonardo, input structure, pseudopotentials, structures from Materials Project |
| 1 | `01_scf_diamond` | first SCF, reading the output |
| 2 | `02_convergence_al` | cutoff, k-point and smearing convergence (fcc Al) |
| 3 | `03_lattice_parameter` | a0 from E(V) and from vc-relax; Pulay stress |
| 4 | `04_lattice_gallery` | how to write fcc, bcc, hcp, diamond, rocksalt, rutile, corundum, wurtzite, layered crystals |
| 5 | `05_fe2o3_afm` | special case: antiferromagnetic oxide, PBE and PBE+U |
| 6 | `06_slab_al111` | slab construction, relaxation, surface energy, orthorhombic cells |
| 7 | `07_adsorption_h_al111` | adsorption energy |
| 8 | `08_neb_h_al111` | energy barriers with NEB: commented template and a worked example |

Also here: `troubleshooting.md` (error messages → causes → fixes),
`tools/` (the Python scripts used by the examples, see `tools/README.md`),
`jobscripts/` (commented Slurm scripts, copied as `job.sh` into each
example).

## Setup (once)

1. Make your working copy (`../README.md`, "Setup").
2. In the copy, put your project budget (see `saldo -b`) in place of
   `CHANGE_ME` in the jobscripts. First list the files, then replace:
   ```bash
   cd qe_runs
   grep -rl CHANGE_ME .
   grep -rl CHANGE_ME . | xargs sed -i 's/CHANGE_ME/<your_project>/'
   ```
   The inputs already point to the group pseudopotential folder
   (`/leonardo/pub/userexternal/epedrett/pseudo`).
3. Load a Python environment for the scripts in `tools/`
   (`../leonardo/README.md`), add the aliases of
   `00_basics/01_running_qe.md` ("Python tools") to your `~/.bashrc`, and
   source it.

## Three habits

1. **Look at every new structure before running it**: `qe2xyz file.pwi`,
   then open the `.extxyz` in OVITO.
2. **Set `max_seconds`** ~5% below the walltime and, for relax and md,
   **`nstep`**: the default is only 50.
3. **Choose `disk_io`**: `'none'` for scf runs whose wavefunctions you do
   not need, `'low'` for relaxations; delete the run's `WORK/` folder (the `outdir`
   of the inputs, not `$WORK`) when a run is finished and you will not
   restart it.

<!-- ## Conventions

- Inputs `.pwi`, outputs with the same name `.pwo`.
- Multi-step calculations: numbered subfolders; the README says what each
  step needs from the previous one.
- Slabs and z-oriented systems: atoms sorted along z; all atoms relaxed
  unless stated.
- Cards in the inputs: `ATOMIC_SPECIES`, `K_POINTS`, `CELL_PARAMETERS`,
  then `ATOMIC_POSITIONS` at the end of the file.
- Reference results were computed with a serial CPU build of QE 7.3.1 and
  the same pseudopotentials; on Leonardo (GPU) expect differences of
  ~1e-6 Ry. -->
