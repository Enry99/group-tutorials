# From Materials Project to a QE input

[Materials Project](https://next-gen.materialsproject.org) (MP) has relaxed
structures for most inorganic crystals, computed with PBE or PBE+U (VASP).
Use it for the starting structure and for a first guess of the k-grid.

## 1. Find the material

Search by formula or elements (e.g. `Fe2O3`). Several entries can exist
for one formula: different polymorphs. Pick by:

- **space group** (α-Fe₂O₃ = corundum, R-3c, No. 167);
- **energy above hull** (0 eV/atom = thermodynamically stable in MP);
- **experimentally observed** (ICSD tag).

Each entry has an ID `mp-NNNN`. Write it in your notes and in the input
header: it is how others find your starting structure.

## 2. Primitive or conventional cell?

The material page lets you export the structure as CIF, POSCAR or JSON,
choosing a **primitive** or **conventional** cell.

| Cell | Atoms (examples) | Use it for |
|---|---|---|
| primitive | bcc Fe: 1, fcc Ni: 1, Fe₂O₃ (rhombohedral): 10 | bulk properties: cheapest |
| conventional | bcc Fe: 2, fcc Ni: 4, Fe₂O₃ (hexagonal): 30 | supercells, slabs along the conventional axes, orthogonal cells |

Same material, same energy per atom (with an equivalent k-point spacing),
different cost. `tools/struct_to_qe.py` prints which one you have and
converts between them.

## 3. Convert to a pw.x input

First see what you downloaded:

```bash
python tools/struct_to_qe.py Fe2O3_mp-XXXX.cif -o Fe2O3.pwi
```

It prints the space group and whether the cell is primitive, and writes
`Fe2O3.pwi` with `ibrav = 0`, `CELL_PARAMETERS angstrom`,
`ATOMIC_POSITIONS crystal`, some default pseudopotentials and default
settings. Options:

- `--primitive` or `--conventional`: convert the cell (spglib);
- `--kspacing 0.20` (default): k-grid from this spacing in 1/Å, a starting
  point to converge (02); `--kgrid 6 6 2` to give it by hand.

Then edit by hand:

- magnetism (`nspin`, `starting_magnetization`, AFM species: see 05);
- DFT+U if the project uses it (MP's U values are for VASP: not transferable);
- the calculation type and thresholds.

Look at it before running: `qe2xyz Fe2O3.pwi` and open the `.extxyz` in
OVITO. Without the script, ASE does the same conversion:
`ase.io.read("file.cif")`, then `ase.io.write("x.pwi", atoms, format="espresso-in", ...)`.

If you prefer `ibrav ≠ 0` (e.g. for `cell_dofree = 'ibrav'`), take the
lattice parameters from MP and the positions from the gallery (04) or
from the CIF, keeping in mind that QE's lattice vectors for each `ibrav`
may be oriented differently from MP's (positions in crystal units refer to
QE's vectors).

## 4. A first guess of the k-grid from the MP calculation

The material ID is also the ID of the original MP calculation (the task
with the numerically smallest ID becomes the material ID). On the material
page, open the list of calculations (*Calculations* / calculation summary)
and find the **GGA (or GGA+U) Structure Optimization** whose task ID is the
`mp-NNNN` of the material: its input shows the KPOINTS used.

How to use those numbers:

- Use MP's grid as a lower bound or initial guess, not as the converged
  value.
- The KPOINTS file says `Gamma` or `Monkhorst-Pack`. In QE use a
  Γ-centered grid (`N N N 0 0 0`) in both cases; a VASP Monkhorst-Pack
  grid with even N is shifted, and shifted grids cause trouble in QE on
  hexagonal cells (`02_convergence_al`, section 5).
- The grid refers to MP's cell. If your cell is different (conventional
  vs primitive, supercell, slab), convert through the k-spacing:
  ```bash
  python tools/kgrid.py POSCAR_from_MP --grid 8 8 8      # spacing of MP's grid
  python tools/kgrid.py my_cell.pwo --kspacing 0.20      # same spacing in my cell
  python tools/kgrid.py POSCAR --kppra 1000              # MP-like density
  ```
- Then converge it (02).

## 5. Compare with MP

- Lattice parameters: MP is PBE (or PBE+U), no D3, VASP PAW. Differences
  of ~0.5-1% with our PBE-D3 results are normal; larger ones mean a problem
  (wrong magnetic state, unconverged cutoff, wrong cell).
- Magnetic ordering: MP lists the ordering it found (FM, AFM, ...) and the
  final moments. A starting point, not a guarantee: for AFM oxides set the
  order explicitly (05).
- Energies: never compare absolute energies across codes or pseudos; only
  differences computed with one set of settings.

Programmatic access (optional): the `mp-api` Python client with a personal
API key from the MP dashboard can download structures and task inputs; see
the MP API documentation. Not needed for the tutorial.
