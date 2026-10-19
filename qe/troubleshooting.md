# Troubleshooting Quantum Espresso

**First rule: look at the structure.** Most errors with a new structure
(SCF that never converges, `S matrix not positive definite`, exploding
forces) are a wrong structure: overlapping atoms, wrong units, wrong
cell. `qe2xyz file.pwi` and open the `.extxyz` in OVITO before anything
else (`00_basics/01_running_qe.md`, "Looking at structures").

Search this page for your error. Errors from pw.x are in the `.pwo`
and in the file `CRASH`; errors from Slurm or MPI are in the job's `.err`.
(Messages marked ✔ were reproduced with pw.x 7.3.1 while writing this page).

New problem? Open an issue in the group repository with the template at
the end, and once solved add it here.

---

## Input errors (pw.x stops at once)

### `file /.../X.UPF not found` (routine readpp) ✔
`pseudo_dir` is wrong, or not set (pw.x then looks in `$HOME/espresso/pseudo/`),
or the file name in `ATOMIC_SPECIES` has a typo. `ls` the path written in
the message.

### `bad line in namelist &system: "..." (error could be in the previous line)` ✔
A misspelled variable, a variable in the wrong namelist (e.g. `ecutwfc`
in `&CONTROL`), or a missing comma/quote on the previous line. Check the
spelling in INPUT_PW.

### `bad line in namelist &ions:` (or `&cell:`) followed by a line of a card ✔
`calculation = 'relax'`/`'vc-relax'` without the `&IONS` (or `&CELL`)
namelist. Add an empty one: `&IONS` + `/`.

### `cannot read atomic specie from: ! ...` / `Error while parsing atomic position card.` ✔
A comment line inside a card. Comments in cards only at the end of a data
line, or between cards.

### `wrong number of columns in ATOMIC_POSITIONS` ✔
`nat` larger than the number of position lines (pw.x read the next card
as a position), or a malformed line.

### `some starting_magnetization MUST be set` ✔
`nspin = 2` with smearing but no `starting_magnetization`. Set it for the
magnetic species.

### `fixed occupations and lsda need tot_magnetization` ✔
`nspin = 2` with `occupations = 'fixed'`. Use smearing, or give
`tot_magnetization`.

### `do not specify both celldm and a,b,c!` ✔
Lattice given both in bohr (`celldm`) and in Å (`A`, `B`, `C`). Keep one. Do not mix parameters from the two sets: e.g. for a hexagonal lattice, use either `celldm(1)` + `celldm(3)` or `A` + `C`, not `celldm(1)` + `C` or `A` + `celldm(3)`.

### `DFT+Hubbard input syntax has changed since v7.1` ✔
Old `lda_plus_u`/`Hubbard_U(i)` variables. Use the `HUBBARD` card (05).

### `ASE does not support ibrav != 0` (Python) ✔
`ase.io.read` on a `.pwi` with `ibrav ≠ 0`. Use `tools/qe2xyz.py` (any
ibrav), or read the `.pwo`.

---
## Memory errors

### The job dies at start-up on the GPU or after a few SCF steps (memory error in `.err`)
Search for `Estimated max dynamical RAM per process` in the `.pwo`. It shows a lower bound. To stay safe, the memory per MPI rank (per GPU) should be at least 1.5× (better 2×) that value. Increase the number of nodes or GPUs, or decrease the k-point parallelization (try to halve `-nk`).

#### Always prefer `K_POINTS gamma` for Γ-point calculations, instead of `automatic` with a 1-point grid (`1 1 1  0 0 0`).
The gamma option uses ~half the memory and is much faster. Exception: DFT+U
with `ortho-atomic` (or `norm-atomic`) projectors, not implemented with
`K_POINTS gamma` (`Gamma-only calculation for this case not implemented`).
With `K_POINTS gamma` and `nspin = 1` use `-nk 1` (`Gamma-only calculations
not allowed with pools`); with `nspin = 2`, `-nk 2` splits the two spins.

---

## SCF problems

### `convergence NOT achieved after N iterations: stopping` ✔
The run ends with `JOB DONE.` anyway, but no energy and forces. In order of what to try:
1. Look at the structure: overlapping or missing atoms, wrong combination of cell and positions (units `crystal`, `alat`, `angstrom`). Atoms too far apart in a bulk indicate a wrong cell. Check the `.pwo` for the last geometry read (or use `qe2xyz`).
2. Decrease `mixing_beta` to 0.1 (or 0.05), `mixing_ndim = 16`, increase `electron_maxstep = 300` (the default is 100 but is often not enough for complex structures).
3. Slabs, interfaces, molecules on surfaces: `mixing_mode = 'local-TF'`.
4. Magnetic systems: check `starting_magnetization` (sign, species); a bad
   start can oscillate between magnetic states.
5. Check the cutoff against the `Suggested minimum cutoff` in the UPF
   header: a pseudopotential used far below it can make the SCF stall
   (seen while testing 05: stuck around 1e-2 Ry with a hard Fe pseudo at
   30 Ry; 40-60 iterations and convergence with a suitable cutoff).
6. As a last resort, just for relaxations, try to set `scf_must_converge = .false.` to attempt overcoming the specific bad atomic structure. Check that all the following steps are fully converged, and inspect the relax trajectory to make sure there are no sudden jumps. Alternatively, first try without spin polarization, then restart with `nspin = 2` and a suitable `starting_magnetization`.

### Isolated atom or molecule in a box: the SCF stalls around 1e-6 Ry ✔
The estimated accuracy oscillates instead of decreasing: the density in
the nearly empty box is numerical noise for the functional (seen with a
spin-polarized H atom in a 12 Å box). For one atom or a small molecule
`conv_thr = 1.0d-6` is enough (0.01 meV); for an isolated atom use fixed
occupations and `tot_magnetization` (07, `01_gas_phase/H_atom.pwi`).

### `S matrix not positive definite` (cdiaghg) / `problems computing cholesky decomposition`
Usually atoms too close (bad structure, interpolation through an atom in
NEB, an MD frame that exploded), or a cutoff too low for the
pseudopotentials. Check the geometry first. If the structure is fine, try
`diagonalization = 'cg'` or `'ppcg'`.

### `too many bands are not converged` (c_bands)
The iterative diagonalization failed for several bands, often in the first
SCF steps of difficult systems. Try `diagonalization = 'cg'` (slower, more
robust); check the structure and magnetization as above.

### The run is magnetic but the result is not (total magnetization 0 for Fe)
`starting_magnetization` missing or too small, or the SCF fell into the
non-magnetic solution. Restart with a larger initial moment.

### AFM system with total magnetization ≠ 0
Wrong spin assignment (species), or the SCF found another magnetic state.
See 05 ("Check the magnetic state").

### A metal looks like an insulator, or `highest occupied level` instead of `the Fermi energy is` ✔
`occupations = 'fixed'` for a metal: pw.x runs without error but the
result is wrong. Use `occupations = 'smearing'`.

### Do not use `smearing = 'methfessel-paxton'` for insulators
It can give wrong occupations (see [VASP wiki](https://vasp.at/wiki/Smearing_technique) for explanation). When in doubt whether the system is metallic (conductive) or insulating, always use Gaussian smearing. Note that "metallic" means that it does not have a gap, not that it contains metal elements.

---

## Relaxations and NEB

### `The maximum number of steps has been reached.` ✔
`nstep` too small for this relaxation. **If `nstep` is not in the input,
the default is only 50** for relax, vc-relax and md: always write it
(e.g. `nstep = 200`). Restart from the last geometry (see "Restart") with
a larger `nstep`.

### `Message from routine bfgs: history already reset at previous step: stopping`
BFGS could not find a lower energy along its direction: usually forces
noisy compared with the threshold. Tighten `conv_thr` (e.g. 1e-8 - 1e-10),
check that `forc_conv_thr` is not stricter than what the SCF accuracy
allows, then restart from the last geometry. If the forces are already
small, the structure may simply be converged within the noise.

### vc-relax: a0 differs from the E(V) fit, or the final scf has a large pressure ✔
Pulay stress: the cutoff is not converged for the stress. The energy
converges much faster: bcc Fe with a hard pseudo at lower cutoff than suggested has a0 from E(V)
within 0.006 Å, while the vc-relax shrinks the cell by 7%, loses the
magnetic moment and ends with −2930 kbar in the final scf. Raise the cutoff
until the two methods agree (03, `03_pulay_fe`).

### Atoms fly away during a relaxation
Forces from an unconverged SCF, or a wrong structure (atoms too close,
wrong cell). A small rigid drift of a slab with all atoms free is
harmless (the net force is zero); large displacements are not.

### NEB: `initial path length` much larger than expected
The two end points do not have the atoms in the same order (e.g. one was
sorted or rebuilt separately). `tools/make_neb.py` stops if the element
sequence differs and warns about large displacements. See 08.

### `Error in routine d_matrix (2): D_S (l=3) for this symmetry operation is not orthogonal` ✔
Met in 08: neb.x stopped at the start of an intermediate image. The
hexagonal cell was not exactly hexagonal: it had been read back from a
`.pwo`, which prints the lattice vectors with 6 digits (relative error
~5e-7). With the cell copied at full precision from the `.pwi` the same
NEB ran. Never build new inputs from the cell of a `.pwo` of an scf or
relax: take it from the input (the scripts in `tools/` do this when the
`.pwi` with the same name is there), or use `ibrav`.

### neb.x with `restart_mode = 'restart'` starts from scratch ✔
No `prefix.path` file in the folder (wrong folder or `prefix`): neb.x
does not stop, it prints `restart_mode = from_scratch` and builds a new
path. Check that line at the top of `neb.pwo`.

---

## Messages that are not errors

### `ecutrho < 4*ecutwfc, are you sure?` ✔
Only a message: pw.x continues. With US/PAW pseudopotentials `ecutrho`
should be 8-12 × `ecutwfc`; check the input.

### `This is a supercell, fractional translations are disabled` ✔
Not an error: the cell is a multiple of a smaller one (e.g. corundum in the
hexagonal cell). pw.x uses fewer symmetry operations, hence more k-points.

### `number of k points=` larger than the grid ✔
A hexagonal cell (slab of fcc(111) or hcp(0001), graphite, MoS₂, MXene)
with a **shifted** grid (`7 7 1 1 1 0`) or with nk1 ≠ nk2, and less
symmetry than the lattice: pw.x adds the points of the grid rotated by
±120°. Example: 7×7×1 = 49 points, 25 without symmetry, but 75 printed.
Cost ×3, and the sampling changes with symmetry and `nosym`, so energies
of different structures do not compare. Use Γ-centered grids with
nk1 = nk2 (`7 7 1 0 0 0`). Details and test: `02_convergence_al`,
section 5.

### `Presently no symmetry can be used with electric field` ✔
Not an error: with the dipole correction (`tefield`/`dipfield`) pw.x uses
no symmetry, so the number of k-points grows.

---

## Time limits and restart

### `Maximum CPU time exceeded` ✔
`max_seconds` reached: pw.x stopped cleanly. To continue:
- `restart_mode = 'restart'` in `&CONTROL`, same input, same `outdir` and
  `prefix` (needs the files in `outdir`, i.e. not `disk_io = 'none'`;
  `'low'` is enough for a relax). Use `>>` instead of `>` before the
  `.pwo` name in the jobscript, otherwise the first part of the output is
  overwritten (`pwo_summary.py` reports the status of the last run in the
  file). The new run prints `Atomic positions and unit cell read from
  directory`; or
- if the restart fails, copy the last geometry from the `.pwo` to a new
  input and start again.

### The job was killed by Slurm (no `JOB DONE.`, `CANCELLED ... DUE TO TIME LIMIT` in `.err`)
`max_seconds` missing or larger than the walltime: no clean stop, the
restart files may be incomplete. Always set `max_seconds` ~5% below the
walltime. Start again from the last positions in the `.pwo`.

---

<!-- ## Template for a new issue

```
Title: <exact error message, first line>

Code/version: pw.x 7.3.1 (module ...), Leonardo booster, N nodes, -nk ...
System: <what>, nat = ..., k-points ..., nspin ...
What I did: <calculation type, what changed from a working input>
Error: <copy the exact lines from .pwo / CRASH / .err>
Files: <path on Leonardo to the run folder>

Solution (when found): <cause> -> <fix>
``` -->
