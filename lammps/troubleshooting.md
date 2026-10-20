# Troubleshooting (LAMMPS)

First the three problems that cost most runs with MLIPs: lost atoms, GPU
out of memory, jobs that stop without a message. Then the other messages,
met in specific inputs: search this page for the exact text. The messages
quoted were reproduced in the local checks (LAMMPS 22 Jul 2025, Update 4;
the group build, 30 Mar 2026, behaves the same).

## Main problems

### `ERROR: Lost atoms: original N current M`

Atoms left the box or moved more than a neighbour-list skin in one step:
the forces were huge. Checklist:

1. Overlapping atoms in the starting structure (stacked parts too close,
   molecules on top of each other): look at it in OVITO, minimize first
   (`05_sliding/01_load.in` does). Reproduced with a molecule film placed
   2.5 Å from a slab and a repulsive potential, at step 0.
2. Timestep too large for H or for the temperature: test it
   (`09_timestep`); 0.5 fs is the usual choice with H.
3. A group integrated twice (warning "One or more atoms are time
   integrated more than once", below), or a driven group with both
   `fix move` and `nve`.
4. Load: force per atom × atoms = P × area?
   (`structures/tools/load_force.py`). A load applied at once on a gap
   makes the top block hit the bottom at high speed: ramp it (`t_ramp`)
   or limit the speed of the approach (below, out of memory).
5. After a restart: the same fixes, in the same order, with their
   `fix_modify` lines (a restart file stores no fixes and no computes;
   give `pair_style` and `pair_coeff` again for the MLIP).
6. A non-periodic boundary (`boundary p p f`): an atom that crosses that
   face is lost. Close it with `fix wall/reflect` (met in the deposition
   template with a wall only at the top).

### GPU out of memory (MLIP)

The cost and the GPU memory of a MACE model grow with the number of atom
pairs within its cutoff (the edges of its graph), not with the number of
atoms. With domain decomposition each GPU holds the pairs of its
subdomain, so what matters is the densest subdomain. Two stages are at
risk:

- **compression / load ramp**: the density and the pairs per atom grow;
- **start of sliding, or a system not yet compact**: atoms cross from one
  subdomain to the next quickly, and one GPU can collect many more pairs
  than at the start.

The run stops with an error containing "out of memory" in the Slurm
output (`grep -i "out of memory" *.out *.err`), or with no message at all
(next section). What helps:

1. **`processors * * 1`** for slabs and interfaces: no cut along z, each
   GPU holds full-height columns, so compression along z does not move
   atoms between subdomains. Balance then in x and y only
   (`balance 1.0 shift xy 5 1.1`, and `fix balance ... shift xy 5 1.1
   weight neigh 0.6` to follow the pairs as they move).
2. **Limit the speed in the compression with the load on**: a longer
   load ramp, or `fix nve/limit` on the moving parts during the approach
   only (it caps the displacement per step; not in production: the
   dynamics is no longer Newtonian).
3. **Headroom**: fewer atoms per GPU (more GPUs) than the first steps
   would allow. Check the peak with a line in the jobscript before
   `mpirun`:
   `nvidia-smi --query-gpu=memory.used,memory.total --format=csv -l 30 > gpu_mem.csv &`
4. **Fragmentation**: the size of the graph changes at every neighbour
   rebuild. `export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` in
   the jobscript lets the PyTorch allocator grow one segment instead of
   leaving free slivers (PyTorch documentation, "CUDA semantics", an
   experimental option). It does not help if the memory is really
   insufficient.

### The job stops without an error

The log ends in the middle of a run, without `ERROR` and without `Total
wall time`. Seen with MACE (ML-IAP), most likely an error inside the
Python/PyTorch part that does not reach the LAMMPS output, probably
memory in most cases.
Find the cause:

1. `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS`:
   `TIMEOUT` (walltime), `OUT_OF_MEMORY` (host memory of the job),
   `NODE_FAIL`, `FAILED` with an exit code.
2. The Slurm `.err` and `.out` files: Python tracebacks and CUDA errors
   end up there, not in the LAMMPS log.
3. The last thermo lines (the templates use `thermo_modify flush yes`):
   at which stage, density and temperature it stopped. A stop during
   compression or at the start of sliding points to memory (section
   above).

Restart from the last restart file (`restart N name.*.restart` in the
templates) with the fixes of the stage, after the change.

## Other messages

### `ERROR: Must use variable energy with fix addforce`

A `minimize` with `fix addforce` driven by a variable and `fix_modify ...
energy yes`. Give the energy of the force as an atom variable:
`variable ez atom -v_fz*z` and `fix ... addforce 0 0 v_fz energy v_ez`.

### `ERROR: Cannot reset timestep with active dump - must undump first`

`reset_timestep` after a `dump` command. Put `reset_timestep` before the
dump (or `undump`, reset, `dump` again).

### `ERROR: Thermo keyword pxx in variable requires thermo to use/init press`

A variable uses `pxx` (or `press`) and `thermo_style` does not print any
pressure. Add `press` to `thermo_style`.

### `ERROR: Variable vx: Invalid special function min() in variable formula`

LAMMPS variables have no `min()` or `max()` of two numbers. Use
`ternary(condition,a,b)`, e.g. `ternary(step<n,step/n,1.0)`.

### `ERROR on proc 0: Potential file ffield... requires real units but metal units are in use`

ReaxFF force fields are in `units real` (kcal/mol, fs). The templates use
`units metal`: times, timestep and energies change with the units.

### `ERROR: Cannot redefine variable ... with a different style`

`-var name value` on the command line defines an index variable; the
input then defines the same name as `equal` or `string`. Only index
variables can be set from the command line (`09_timestep`).

### `WARNING: One or more particle depositions were unsuccessful`

`fix deposit` found no place in the insertion region at least `near` Å
from the other atoms, usually because the growing film has reached the
region. Fewer atoms are inserted than asked. Move the region higher (more
vacuum in the data file) or deposit fewer atoms (`08_deposition`).

### `WARNING: One or more atoms are time integrated more than once`

Two integrators (`nve`, `nvt`, `move`, `rigid`, ...) act on the same
atoms: they move twice per step. Treat it as an error. Typical: `fix nve`
on `all` plus `fix move` on a layer; `nve` on a union group that contains
a thermostatted group; `fix rigid` and `nve` on the same layer.

### `WARNING: Temperature for fix modify is not for group all`

Printed for every thermostat that uses its own `temp/com` compute
(`fix_modify ... temp`), as in the sliding templates. Expected.

### `WARNING: Cannot count rigid body degrees-of-freedom before bodies are fully initialized`

Once, at the start of a run with `fix rigid`. Does not affect the run.

### `WARNING: Energy due to 1 extra global DOFs will be included in minimizer energies`

`minimize` with `fix box/relax`: the cell is a degree of freedom. Expected.

### Silent problems

- **Slab cut in two in the unwrapped trajectory** after opening a vacuum
  gap in a bulk from MD: image flags from the bulk run. `set group all
  image 0 0 0` after `change_box` (`03_meltquench/anneal_slab.in`;
  `image NULL NULL 0` is ignored by 22 Jul 2025 up to Update 2).
- **A loop that never ends** (`label` ... `jump SELF label` without `next`):
  the run goes on until the walltime, the lines after the loop never run.
  Found in older group CVD inputs.
- **`fix move` units**: the default is `units lattice`. Write `units box`
  (with no `lattice` command the two coincide, but the next input may
  have one).
- **`fix_modify` lost**: a fix redefined with the same ID replaces the old
  one, and its `fix_modify` settings are gone: repeat them.
