# 08 · Deposition

Atoms or molecules are inserted above a surface at regular intervals and
shot at it, as in CVD or PVD growth; the film grows atom by atom. In
cycles, the atoms that bounced back or desorbed (free atoms high above the
surface) are deleted, so that they do not accumulate under the top wall.
You get the grown film and the sticking statistics (inserted vs deleted).
From the group LPCVD/PECVD runs of SiNₓ (SiH₄ and NH₃, removal of H₂ and
N₂).

The test system: `Si001.data`, Si(001) with 8 layers (128 atoms) and 30 Å
of vacuum, groups in `Si001.groups.lmp`. Built with:

```bash
python -c "from ase.build import diamond100; from ase.io import write; write('Si001.extxyz', diamond100('Si', size=(4,4,8), a=5.431, vacuum=10))"
python ../../structures/tools/interface.py Si001.extxyz --no-top --fixed 2 --thermo 4 --vacuum 30 --specorder Si -o Si001.data
```

```bash
sbatch job.sh deposition.in
```

- `fix deposit`: `ncycle × nper` insertions of type `type_dep`, one every
  `every` steps, in a 5 Å slab `zin` Å above the surface, with velocity
  `−vz` along z. The thermal speed √(3k_BT/m) is 9.4 Å/ps for Si at 1000 K;
  higher speeds model energetic (PVD, ion) deposition.
- Molecules: `molecule` command with a molecule file and `mol` in `fix
  deposit` (commented lines in the input; molecule file format in the
  LAMMPS manual, command `molecule`).
- Every cycle (`nper` insertions): `delete_atoms group free`, where `free`
  is a dynamic group of atoms with no neighbour within `rcut` and more than
  5 Å above the initial surface. For products such as H₂, use the
  coordination with the commented lines (`07_reactive` counts them the
  same way).
- The loop: `variable icycle loop`, `label`, `run`, `next icycle`, `jump
  SELF cycle`. The `next` line is what ends it: without it the loop runs
  until the walltime and the final `write_data` is never reached.

The `atoms` column of the thermo output gives the number of atoms; the
log reports `Deleted N atoms` at each cycle.

## Local check (reduced settings)

`pair_style tersoff` with `Si.tersoff` (LAMMPS `potentials/`), 5 cycles of
10 insertions (100 ps), serial CPU build, 21 s: 50 atoms inserted, 1 atom
deleted per cycle (bounced back), 173 atoms at the end; thermostatted
layer between 830 and 1220 K (128 atoms: large fluctuations). Log in `reference/`.
