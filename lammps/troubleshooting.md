# Troubleshooting (LAMMPS)

Search this page for the exact text of the message. ✔: reproduced with
LAMMPS 22 Jul 2025 (Update 4), the version of the group build.

## Errors

### `ERROR: Lost atoms: original N current M` ✔

Atoms left the box or moved more than a neighbour-list skin in one step:
the forces were huge. Checklist:

1. Overlapping atoms in the starting structure (stacked parts too close,
   molecules on top of each other): look at it in OVITO, minimize first
   (`05_sliding/01_load.in` does). Reproduced with a molecule film placed
   2.5 Å from a slab and a repulsive potential, at step 0.
2. Timestep too large: 0.5 fs with H or at high T.
3. A group integrated twice (see the warning below), or a driven group
   with both `fix move` and `nve`.
4. Load: force per atom × atoms = P × area? (`structures/tools/load_force.py`)
5. After a restart: the same fixes, in the same order, with their
   `fix_modify` lines.

### `ERROR: Must use variable energy with fix addforce` ✔

A `minimize` with `fix addforce` driven by a variable and `fix_modify ...
energy yes`. Give the energy of the force as an atom variable:
`variable ez atom -v_fz*z` and `fix ... addforce 0 0 v_fz energy v_ez`.

### `ERROR: Cannot reset timestep with active dump - must undump first` ✔

`reset_timestep` after a `dump` command. Put `reset_timestep` before the
dump (or `undump`, reset, `dump` again).

### `ERROR: Thermo keyword pxx in variable requires thermo to use/init press` ✔

A variable uses `pxx` (or `press`) and `thermo_style` does not print any
pressure. Add `press` to `thermo_style`.

### `ERROR: Variable vx: Invalid special function min() in variable formula` ✔

LAMMPS variables have no `min()` or `max()` of two numbers. Use
`ternary(condition,a,b)`, e.g. `ternary(step<n,step/n,1.0)`.

### `ERROR on proc 0: Potential file ffield... requires real units but metal units are in use` ✔

ReaxFF force fields are in `units real` (kcal/mol, fs). The templates use
`units metal`: times, timestep and energies change with the units.

## Warnings

### `WARNING: One or more atoms are time integrated more than once` ✔

Two integrators (`nve`, `nvt`, `move`, `rigid`, ...) act on the same
atoms: they move twice per step. Treat it as an error. Typical: `fix nve`
on `all` plus `fix move` on a layer; or `nve` on a union group that
contains a thermostatted group.

### `WARNING: Temperature for fix modify is not for group all` ✔

Printed for every thermostat that uses its own `temp/com` compute
(`fix_modify ... temp`), as in the sliding templates. Expected.

### `WARNING: Cannot count rigid body degrees-of-freedom before bodies are fully initialized` ✔

Once, at the start of a run with `fix rigid`. Does not affect the run.

### `WARNING: Energy due to 1 extra global DOFs will be included in minimizer energies` ✔

`minimize` with `fix box/relax`: the cell is a degree of freedom. Expected.

## Silent problems

- **Slab cut in two in the unwrapped trajectory** after opening a vacuum
  gap in a bulk from MD: image flags along z from the bulk run.
  `set group all image NULL NULL 0` after `change_box`
  (`03_meltquench/anneal_slab.in`). ✔
- **A loop that never ends** (`label` ... `jump SELF label` without `next`):
  the run goes on until the walltime, the lines after the loop never run.
  Found in older group CVD inputs.
- **`fix move` units**: the default is `units lattice`. Write `units box`
  (with no `lattice` command the two coincide, but the next input may
  have one).
