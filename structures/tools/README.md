# tools/

Scripts used by the examples. They need ase, numpy and, where noted,
pymatgen or packmol (`../README.md`, "Setup"); `-h` prints the usage of each
one. All of them write `.extxyz` by default, other formats with `-f`
(`-f lammps-data --specorder ...` for LAMMPS).

| Script | What it does |
|---|---|
| `slab.py` | slabs of any bulk and Miller index, one file per termination; stoichiometry and top/bottom planes printed; c along z, vacuum, supercell (pymatgen) |
| `zsl_interface.py` | coherent interface of two crystals by lattice matching: list of matches (area, strain), chosen interface and the two strained slabs (pymatgen) |
| `density.py` | repeat a bulk and scale it to a target density; `--xy` keeps the lateral cell of the slab it will sit on |
| `dope.py` | random substitutions (fraction of all atoms) and interstitials (covalent-radius check) |
| `packmol.py` | molecules in a box (or the lateral cell of a slab) at a target density, optionally around a fixed solute (packmol) |
| `interface.py` | bottom slab / middle / top slab or tip → LAMMPS data, `<name>.groups.lmp` with the groups, force per atom for a load |
| `tip.py` | hemispherical tip from a bulk, optional H/OH passivation of the curved surface |
| `load_force.py` | pressure (GPa) ↔ force per atom (eV/Å) for `fix addforce` |
| `common.py` | helpers shared by the scripts (sorting along z, density, writing LAMMPS data) |

Also used: `../../qe/tools/ortho_cell.py`, oblique in-plane cell →
smallest orthorhombic supercell (`qe/06_slab_al111`, "Orthorhombic cells").
