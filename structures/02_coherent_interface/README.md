# 02 · Coherent interfaces (Zur-McGill lattice matching)

An interface between two different crystals is periodic only if their
surface lattices fit: supercells of the two that coincide within a small
strain. The Zur-McGill algorithm (A. Zur and T. C. McGill, J. Appl. Phys.
55, 378 (1984), doi:10.1063/1.333084) lists them; pymatgen implements it
(`ZSLGenerator`) and builds the interface (`CoherentInterfaceBuilder`). Use it
for film/substrate and heterogeneous contacts (oxide on metal, nitride on
diamond). `tools/zsl_interface.py` lists the matches and writes the one you
choose. You get the interface and the two strained slabs.

Example: rutile TiO₂(110) on bcc Fe(110). The film (TiO₂) is strained onto
the substrate (Fe).

## 1. List the matches

```bash
python ../tools/zsl_interface.py ../bulk/TiO2_rutile.cif 1 1 0 ../bulk/Fe_bcc.cif 1 1 0 --list
```

```
terminations (film, substrate):
  0: ('TiO_Pmmm_4', 'Fe_Cmmm_1')
  1: ('O2_Pmmm_1', 'Fe_Cmmm_1')
termination 0: 389 matches, 97 distinct
    i   a (A)    b (A)  gamma   area (A^2)  strain (%)  atoms
    0    7.020   11.007   77.7        75.5       1.30    111
    ...
```

The terminations are those of `slab.py` (01): termination 1 is the
O-terminated TiO₂(110). Matches are sorted by area. For MD prefer
γ = 90° (sliding along x or y) and low strain:

```bash
python ../tools/zsl_interface.py ../bulk/TiO2_rutile.cif 1 1 0 ../bulk/Fe_bcc.cif 1 1 0 \
       --termination 1 --list | awk 'NR <= 5 || $4 == 90.0'
```

```
    7    7.020   32.266   90.0       226.5       1.30    333
    9   14.330   16.213   90.0       232.3       0.83    336
   11   12.159   20.062   90.0       243.9       1.92    342
```

Tolerances: `--length-tol 0.03` and `--angle-tol 0.01` (fractions),
`--area-tol 0.09`, `--max-area 400` Å². Larger values give more matches
with more strain. The search takes about 15 s here.

## 2. Write the interface

```bash
python ../tools/zsl_interface.py ../bulk/TiO2_rutile.cif 1 1 0 ../bulk/Fe_bcc.cif 1 1 0 \
       --termination 1 --select 9 --film-layers 3 --substrate-layers 3
```

```
match 9: 336 atoms (216 film, 120 substrate), a = 14.3300 b = 16.2125 A, gamma = 90.00 deg, film strain (von Mises) 0.83 %
```

It writes `TiO2_rutile_Fe_bcc_9.extxyz` (Fe at the bottom, 2.5 Å gap
(`--gap`), TiO₂ above, 15 Å vacuum), and the two strained slabs
`..._film.extxyz` and `..._substrate.extxyz`. Checks on the result: Fe-Fe
distance 2.482 Å as in the bulk (the substrate is not strained), shortest
Ti-O 1.946 Å.

Thicknesses are in unit planes: 3 Fe(110) planes are 4.05 Å, too thin for
anything but a test. Use the separate slabs to build a thicker system or a
sliding interface with `interface.py` (04).

## Notes

- The strain goes on the film. Choose as substrate the material whose
  lattice you want to keep (the thicker or stiffer one).
- 0.8 % strain on a film changes its energy and can drive artefacts in MD
  (e.g. misfit dislocations that are only a cell effect). Compare a few
  matches before production.
- The gap is only a starting value: relax the interface (or minimize under
  a small load in LAMMPS) before anything else.
