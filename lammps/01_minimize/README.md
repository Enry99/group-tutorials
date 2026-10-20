# 01 · Energy minimization

Minimization moves the atoms (and, with `box/relax`, the cell) downhill in
energy until the forces are small. Use it before MD on a structure you
built (overlaps, a gap between stacked parts), to get the lattice parameter
of the potential at 0 K, and at the end of a run to get a structure for
analysis or DFT. You get the relaxed structure (`minimized.data`) and, with
box/relax, the 0 K cell.

`minimize.in` (on `Fe_bcc.data`, 250 atoms of bcc Fe, also used by 02 and
04; built with `python -c "from ase.build import bulk; from ase.io import
write; write('Fe_bcc.data', bulk('Fe', 'bcc', a=2.866, cubic=True).repeat(5),
format='lammps-data', masses=True)"`):

1. atoms only, fixed cell: `minimize 0.0 1.0e-3 10000 100000`, i.e. no
   energy criterion, largest force component below 1e-3 eV/Å
   (`min_modify norm inf`), at most 10000 iterations;
2. atoms and cell at zero pressure: `fix box/relax iso 0.0`, then the same
   `minimize`. Use `aniso` for non-cubic cells, `x 0.0 y 0.0` for slabs
   (never relax z across the vacuum).

```bash
sbatch job.sh minimize.in
grep -E "Stopping criterion|final cell" minimize.log    # when the job has finished
```

What to check: `Stopping criterion = force tolerance` (not `max
iterations`), and `fmax` in the last thermo line below the tolerance. With
an MLIP in single precision the minimizer can stop on `linesearch alpha is
zero` before the tolerance: then look at `fmax`; a few 1e-3 eV/Å is
usually as low as it goes. With box/relax the pressure at the end should
be close to the target (a few bar).

## Notes

- `WARNING: Energy due to 1 extra global DOFs will be included in
  minimizer energies` with box/relax: expected (the cell is a degree of
  freedom).
- Forces added by `fix addforce` (a load) enter a minimization only with
  `fix_modify ID energy yes` and, for a variable force, the `energy`
  keyword of addforce (`05_sliding/01_load.in`).
- A minimization does not cross barriers: an amorphous or strained
  structure stays in the nearest minimum. For that use MD (anneal) first.

## Local check (reduced settings)

`pair_style eam/fs` with `Fe_mm.eam.fs` (Mendelev EAM, in the LAMMPS
`potentials/` folder) in place of the MLIP, serial CPU build: the perfect
crystal is already at a force minimum (1 iteration), box/relax gives
a = 14.2766/5 = 2.8553 Å in 4 iterations, P = 0.1 bar.
Log: `reference/minimize_cpu-serial.log`.
