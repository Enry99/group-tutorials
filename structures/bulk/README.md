# bulk/

Bulk cells used by the examples, written with ASE (`ase.spacegroup.crystal`)
from the experimental parameters of `qe/04_lattice_gallery` and
`qe/01_scf_diamond`:

| File | Structure | Parameters |
|---|---|---|
| `C_diamond.cif` | diamond, Fd-3m (227) | a = 3.567 Å |
| `Fe_bcc.cif` | bcc, Im-3m (229) | a = 2.866 Å |
| `TiO2_rutile.cif` | rutile, P4₂/mnm (136) | a = 4.594 Å, c = 2.959 Å, u = 0.3048 |

For production, start from the bulk relaxed with the method you will use
(DFT settings or the MLIP): a slab cut from the experimental cell is
strained for that method.
