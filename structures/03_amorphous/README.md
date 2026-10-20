# 03 · Starting structures for amorphous materials

Amorphous solids (a-C, DLC, a-SiNₓ:H, oxide films) are made in MD by
melt-quench: a crystal at the right composition and density is melted,
equilibrated as a liquid and cooled (`lammps/`, to come). The starting
structure only has to fix composition, number of atoms and density: the
melt forgets the crystal. Here you build it from a bulk crystal by
substitution, insertion and rescaling.

- `tools/dope.py`: random substitutions (fraction of all atoms) and
  interstitials (with a check on covalent radii).
- `tools/density.py`: repetition and scaling to the target density, on all
  sides or only along z (`--xy`, for a film on a given slab).

## 1. Doped DLC: Si and H in carbon at 3.0 g/cm³

```bash
python ../tools/dope.py ../bulk/C_diamond.cif --repeat 5 4 4 --substitute Si:0.10 H:0.20
python ../tools/density.py C_diamond_Si10_H20.extxyz --target 3.0 -f lammps-data --specorder C Si H
```

```
start: C640 (640 atoms), 3.516 g/cm^3
substituted 64 host atoms with Si
substituted 128 host atoms with H
final: C448H128Si64 (640 atoms), 3.342 g/cm^3
...
final: 640 atoms, cell 18.4886 x 14.7909 x 14.7909 A, density 3.000 g/cm^3
-> C_diamond_Si10_H20_3.00gcm3.data
```

The fractions are of all atoms (10 % Si = 64 of 640). `--seed` changes the
random choice: use different seeds for independent samples.

## 2. Interstitial H: lower the density first

Interstitials need room. In diamond at 3.5 g/cm³, 10 % H does not fit
(`only 7 of 51 H placed after 20000 attempts`, about 20 s): scale first,
then insert:

```bash
python ../tools/density.py ../bulk/C_diamond.cif --repeat 4 4 4 --target 2.8 -o C_2.8.extxyz
python ../tools/dope.py C_2.8.extxyz --interstitial H:0.10
```

```
inserted 51 H (2676 attempts)
final: C512H51 (563 atoms), 2.823 g/cm^3
```

A new atom is rejected if it is closer than 1.4 × (sum of covalent radii)
to any atom (`--mult`, the value of `si3n4_doping.py`).

## 3. Film for a given slab: keep the lateral cell

A film that will sit on a slab (or between two slabs) must have the slab's
lateral cell. `--xy` sets it and scales only c to reach the density. Here
TiO₂ for the 19.856 × 21.061 Å Fe(110) slab of `01_slabs`. 3.8 g/cm³ is
only an example value: use the density of the amorphous phase you want,
from experiment or literature.

```bash
python ../tools/density.py ../bulk/TiO2_rutile.cif --repeat 4 3 6 --xy 19.8562 21.0607 \
       --target 3.8 -f lammps-data --specorder Fe Ti O -o TiO2_for_Fe110.data
```

```
start: 432 atoms, cell 18.3760 x 13.7820 x 17.7540 A, density 4.247 g/cm^3
lateral cell set: 432 atoms, cell 19.8562 x 21.0607 x 17.7540 A, density 2.572 g/cm^3
final: 432 atoms, cell 19.8562 x 21.0607 x 12.0175 A, density 3.800 g/cm^3
```

The crystal is now strongly strained (b stretched by 53 %, c compressed):
irrelevant for a melt-quench. `--specorder Fe Ti O` already contains Fe,
so the same type numbering works when the film is put on the Fe slab.

Choose `--repeat` so that the number of atoms gives the cell you need: N
fixes the volume at a given density. The full outputs are in `reference/`.

## Notes

- **Density.** It controls the sp³ fraction of a-C, the coordination of
  a-SiNₓ, the porosity of oxides. Check it after the quench: an NVT quench
  keeps it, an NPT quench lets the potential choose it.
- **Composition** is fixed at this step. H added as substitution keeps the
  number of atoms; as interstitial it adds atoms.
- The melt-quench itself (temperatures, rates, NVT or NPT) is in the LAMMPS
  part.
