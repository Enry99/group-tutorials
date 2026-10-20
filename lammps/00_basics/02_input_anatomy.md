# Anatomy of a LAMMPS input

A LAMMPS input is a script: commands are executed one after the other,
and `run` or `minimize` advance the system. The templates follow the same
order of blocks. Reference for every command:
[LAMMPS manual](https://docs.lammps.org/Commands_all.html).

## Blocks

```
# ---------- settings ----------      variable model, elements, T, times
# ---------- system ----------        units, atom_style, boundary, read_data / read_restart, groups
# ---------- potential ----------     pair_style, pair_coeff
# ---------- integration ----------   timestep, computes, fixes (integrators, thermostats, forces)
# ---------- output ----------        thermo, dump, restart
# ---------- 1. stage ----------      velocity, run, write_restart
```

- `units`, `atom_style`, `boundary` come before `read_data`. After
  `read_restart` they are already set (they are in the restart file).
- **Fixes act in every following `run`** until `unfix`. A fix defined again
  with the same ID replaces the old one, but its `fix_modify` settings are
  lost: repeat them (as in `07_reactive/highT.in`).
- Every atom that moves needs one integrator (`nve`, `nvt`, `npt`,
  `rigid`, `move`), and only one. Atoms with no integrator stay still.

## Units (`units metal`)

| Quantity | Unit | Useful |
|---|---|---|
| distance | Å | |
| time | ps | `timestep 0.001` = 1 fs |
| energy | eV | |
| velocity | Å/ps | 1 m/s = 0.01 Å/ps |
| force | eV/Å | 1 eV/Å = 1.602 nN |
| pressure, stress | bar | 1 GPa = 10⁴ bar |
| temperature | K | |
| density (thermo `density`) | g/cm³ | |

## Variables

```
variable  T      equal  300                    # a number
variable  model  string CHANGE_ME.model-mliap_lammps.pt
variable  ztop   equal  xcm(layer_top,z)      # evaluated when used (v_ztop)
variable  ez     atom   -v_fz*z               # one value per atom
```

- `${T}` is replaced by its value when the line is read; `v_T` is
  evaluated every time it is used (in a fix, at every step).
- `$(expression)` computes on the spot: `run $(round(v_t_hold/dt))`.
- An `equal` variable can be redefined (`variable fz equal ...` again):
  the fixes that use `v_fz` take the new formula from the next step.
- Loops: `variable i loop 10`, `label name`, ..., `next i`,
  `jump SELF name`. Without `next` the loop never ends (`08_deposition`).

## Restart files

- `write_restart file` stores atoms, velocities, cell, groups, and the
  state of some fixes. It does **not** store fixes, computes, variables,
  dumps, and for the MLIP pair styles not even the potential: the input
  that reads it defines them again (`05_sliding/02_slide.in`).
- `restart N name.*.restart` writes one every N steps (`*` = step).
- `write_data file nocoeff`: a data file with positions and velocities,
  without the pair coefficients; readable by `read_data`, OVITO and ASE.
- `reset_timestep 0` restarts the step counter: time, ramps written with
  `step` and analysis time axes start again from 0. Not allowed with an
  active dump: put it before the `dump` command.

## Outputs

- `thermo N` + `thermo_style custom ...`: a line every N steps in the log
  (`.log` with the jobscripts). Groups and fixes give their own columns:
  `c_ID` (compute), `f_ID[i]` (fix), `v_name` (variable).
- `dump ID all custom N file.lammpstrj id element xu yu zu` +
  `dump_modify ID element <elements> sort id`: the trajectory for OVITO.
  `xu yu zu` are unwrapped: an atom that crosses the cell keeps moving
  continuously (needed for diffusion and velocity profiles).
- Profiles along z during the run: `compute chunk/atom bin/1d` +
  `fix ave/chunk` (e.g. density of each species across a sliding gap).

## Before submitting

1. `variable elements` in the order of the atom types of the data file
   (`Masses` section of the data file).
2. `variable model`: the path exists, and the model has all the elements.
3. Timestep: 1 fs; 0.5 fs with H or above ~1500 K. Check energy
   conservation once (`02_md_ensembles/nve.in`).
4. Every moving group has exactly one integrator (no warning `One or more
   atoms are time integrated more than once`).
5. A short test in the debug queue before the long run.
