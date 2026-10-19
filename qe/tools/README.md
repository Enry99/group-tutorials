# tools/

Python scripts used by the examples. All need the Python environment of
`../../leonardo/README.md` (ase, numpy, spglib, matplotlib); `-h` prints
the usage of each one. The inputs they write start with a comment line
`! written by: <command>`, so you can always see how a file was made.

**Every day**

| Script | What it does |
|---|---|
| `qe2xyz.py` | `.pwi` (any ibrav, NEB inputs), `.pwo` (all ionic steps), `neb.axsf` (all NEB images) → `.extxyz` for OVITO, or cif, xsf, lammps-data; keeps QE species, fixed atoms, starting magnetization; warns about overlapping atoms |
| `pwo_summary.py` | status of `.pwo` files (DONE, SCF-FAIL, TIMEOUT, MAXSTEP, CRASHED), energy, forces, pressure, magnetization, wall time |

**Bulk (02-04)**

| Script | What it does |
|---|---|
| `make_scan.py` | series of inputs from a template: cutoff, k-grid, lattice scaling, any namelist variable (`--set` fixes others, e.g. the cutoff of an E(V) scan) |
| `conv_plot.py` | convergence table and plot of a scan |
| `eos_fit.py` | Birch-Murnaghan (or other) E(V) fit with ASE: V0, B0, a0 |
| `kgrid.py` | k-grid ↔ k-spacing; Materials Project k-point density |
| `struct_to_qe.py` | CIF/POSCAR (e.g. from Materials Project) → pw.x input with default settings; primitive/conventional cell |

**Surfaces, adsorption, NEB (06-08)**

| Script | What it does |
|---|---|
| `build_slab.py` | slab from the relaxed bulk (fcc, bcc, hcp surfaces), settings copied from the bulk input, atoms sorted along z |
| `surface_energy.py` | γ from slab and bulk outputs |
| `add_adsorbate.py` | one atom on a site of the relaxed slab (`--offset` for an equivalent site in the next surface cell), settings copied from the slab input |
| `ads_energy.py` | E_ads = E(slab+X) − E(slab) − E(X) from the three outputs; checks cutoffs and atom count |
| `make_neb.py` | neb.x input from two relaxed end points (`--via` for intermediate images); checks the atom order |
| `ortho_cell.py` | oblique in-plane cell (hexagonal, bcc(110), ...) → smallest orthorhombic supercell, in-plane repetition, LAMMPS data |

**Other files**

| File | What it is |
|---|---|
| `group_pseudos.json` | default pseudopotential per element (`struct_to_qe.py`, `add_adsorbate.py`) |
| `qe_common.py` | helpers shared by the scripts (loading structures, sorting along z, reading templates) |
