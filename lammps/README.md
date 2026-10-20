# LAMMPS templates

Commented LAMMPS inputs for the simulations we run with machine-learning
potentials: relaxation, MD in the basic ensembles, melt-quench, deformation,
sliding interfaces, indentation and scratch, reactions at high temperature,
deposition, and a test of the timestep. This is a catalogue rather than a
course: each folder is one type of simulation. Copy the template closest
to what you need, then change the settings block at the top. For the
commands themselves see the [LAMMPS manual](https://docs.lammps.org).

(Group build on Leonardo: LAMMPS 30 Mar 2026, development version
16f8568, with Kokkos and ML-IAP; see `00_basics`.)

## Contents

| # | Folder | Template | What it does |
|---|---|---|---|
| 0 | `00_basics/` | | running on Leonardo, structure of an input, MLIP pair styles, restarts, outputs |
| 1 | `01_minimize` | `minimize.in` | energy minimization, atoms and cell |
| 2 | `02_md_ensembles` | `nve.in`, `nvt.in`, `npt.in` | MD at constant energy, temperature, pressure |
| 3 | `03_meltquench` | `meltquench.in`, `anneal_slab.in` | amorphous bulk from a crystal; its surface |
| 4 | `04_deformation` | `uniaxial.in` | stress-strain curve at constant strain rate |
| 5 | `05_sliding` | `01_load.in`, `02_slide.in` | slab / lubricant / slab: load, then sliding at constant load |
| 6 | `06_indentation` | `01_indent.in`, `02_scratch.in` | tip on a slab: indentation (displacement control), scratch at constant load |
| 7 | `07_reactive` | `highT.in` | molecules on a surface at high T, molecule counts during the run |
| 8 | `08_deposition` | `deposition.in` | atoms or molecules shot at a surface, removal of the free ones |
| 9 | `09_timestep` | `equil.in`, `nve.in`, `layout.in` | timestep test with H: energy conservation, temperatures, H chemistry vs dt; `check_dt.py` |

The data files of the templates are built with `../structures/` (each
README gives the command). Also here: `jobscripts/` (Slurm scripts, copied
as `job.sh` into each folder) and `troubleshooting.md`.

## Setup (once)

1. Make your working copy (`../README.md`, "Setup"):
   `cp -r group-tutorials/lammps lammps_runs`.
2. Put your project budget (`saldo -b`) in place of `CHANGE_ME` in the
   jobscripts and the path of your model in the inputs (`variable model`,
   also `CHANGE_ME`). List the files first:
   ```bash
   cd lammps_runs
   grep -rl CHANGE_ME .
   ```
   The account can be replaced everywhere at once with
   `sed -i 's/--account=CHANGE_ME/--account=<your_project>/' */job.sh jobscripts/*.sh`;
   the model is usually different for each project, set it by hand.
3. Set `variable elements` in each input to the element list of the model,
   in the order of the atom types of your data file
   (`../structures/README.md`, "Conventions").

## Conventions of the templates

- Settings at the top (`variable ...`): model, elements, temperatures,
  times in ps. The number of steps is computed from the times and the
  timestep (`$(round(v_t_hold/dt))`), so changing `timestep` keeps the
  physical times.
- `units metal` (Å, ps, eV, bar, K), `atom_style atomic`.
- Thermostats: Nose-Hoover, `Tdamp = 100 dt`; barostats `Pdamp = 1000 dt`.
  Regions that move (sliding, indentation) are thermostatted with the
  temperature without the centre-of-mass motion (`compute temp/com` and
  `fix_modify ... temp`).
- Trajectories: `dump custom ... id element xu yu zu` with `sort id`:
  element names (OVITO, ASE) and unwrapped coordinates (velocity
  profiles, diffusion).
- One `write_restart` at the end of each stage, so a long run can be
  continued or branched.

## Local checks

The templates were run in a container with a serial CPU build of LAMMPS
22 Jul 2025 (Update 4) and classical potentials in place of the MLIP (EAM
for Fe, Tersoff for C and Si, ReaxFF for Fe/O/H), at reduced settings.
This checks the inputs, not the physics: each README says which potential
and settings were used; logs, plots and pictures are in `reference/`
(`*_cpu-serial.*`). The changes between that version and the group build
do not touch the commands of the templates (LAMMPS git history and manual
pages checked). A bug of older 2025 versions in `set image` is avoided in
`03_meltquench/anneal_slab.in`.
