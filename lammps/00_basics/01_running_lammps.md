# Running LAMMPS on Leonardo

## The group build

```
/leonardo/pub/userexternal/epedrett/software/lammps-22Jul2026/build/lmp
```

LAMMPS 30 Mar 2026, development version (git commit 16f8568), built with
Kokkos for the A100 GPUs and the ML-IAP package with Python, which runs
MACE models (`03_mlip_pair_styles.md`). Check the version and the packages
of any binary with `lmp -h | head -40`. It needs:

```bash
module purge
module load cuda/12.6 openmpi/4.1.6--gcc--12.2.0
source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate   # torch, mace
```

The command line switches Kokkos on:

```bash
mpirun -np 1 $LMP -k on g 1 -sf kk -pk kokkos newton on neigh half -in input.in
```

- `-k on g 1`: Kokkos with 1 GPU per node (`g 4` on a full node);
- `-sf kk`: use the Kokkos version of every style that has one;
- `-pk kokkos newton on neigh half`: the settings ML-IAP needs.

One MPI rank per GPU. With the jobscripts you never type this line.

## The jobscripts

In `../jobscripts/`, copied into each folder as `job.sh`:

| Script | Resources | Used for |
|---|---|---|
| `lmp_1gpu.sh` | 1 GPU, 8 cores, debug queue (30 min, max 2 jobs) | the templates; MLIP systems up to ~10⁴-10⁵ atoms |
| `lmp_node.sh` | 4 GPUs per node, `--nodes` N | large systems (not tested on the templates) |

```bash
sbatch job.sh 01_load.in        # -> 01_load.out (screen) and 01_load.log (LAMMPS log)
```

For production runs delete the `--qos=boost_qos_dbg` line and set
`--time` (max 24 h). Runs longer than 24 h are split into restarts.

## Chained jobs

A run that needs the restart of another one can be queued at once:

```bash
id=$(sbatch --parsable job.sh 01_load.in)
sbatch --dependency=afterok:$id job.sh 02_slide.in
```

`afterok`: the second job starts only if the first one ended without
error. To continue the same stage after the walltime, start from the last
periodic restart (`restart N name.*.restart` in the templates) with an
input that begins with `read_restart`.

## Checks on a run

- The `.log` header: `KOKKOS mode with Kokkos version ... is enabled`,
  and the GPU that is used.
- Speed: `Performance: ... ns/day, ... timesteps/s` at the end of each run.
  Compare it when you change the number of GPUs or nodes: check that more
  GPUs really make the run faster for your system size.
- `WARNING: Temperature for fix modify is not for group all`: expected in
  the templates (a thermostat on a region with its own temperature compute).
- With `-sf kk`, styles without a Kokkos version run on the CPU and copy
  data between CPU and GPU at every step. In the templates: `fix move`,
  `fix rigid`, `fix deposit`, `compute reduce`. With a MLIP the model
  dominates the cost, but check the timing breakdown at the end of the log
  (`Modify` line) if a run is slower than expected.
