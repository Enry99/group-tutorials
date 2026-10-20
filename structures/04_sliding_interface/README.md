# 04 · Sliding interface: slab / lubricant / slab

The standard tribology cell: two slabs, with molecules (lubricant, water,
additives) or a film in between, periodic in x and y. In LAMMPS the outer
layers of each slab hold the system: the lowest layer of the bottom slab
is fixed, the outer layer of the top slab carries the load and the sliding
motion, the next layers are thermostatted, the rest evolves freely
(NVE). `tools/interface.py` stacks the parts, writes the LAMMPS data file
and a file with the group commands, and computes the force per atom for a
given load.

Example: Fe(110) / water with 10 O₂ / Fe(110).

## 1. Slab and molecules

```bash
python ../tools/slab.py ../bulk/Fe_bcc.cif --miller 1 1 0 --layers 6
python ../../qe/tools/ortho_cell.py Fe_bcc_110_t0.extxyz --square 20 -o Fe110.extxyz
python ../tools/packmol.py H2O O2:10 --cell-from Fe110.extxyz --thickness 15 --density 1.0 -o water
```

```
box [19.856 21.061 15.   ] A, free volume 6273 A^3: 191 H2O, 10 O2
packmol: Success! 593 atoms -> water.extxyz
```

The first molecule (H₂O) fills the box at 1.0 g/cm³ together with the
others (10 O₂). Molecules can also be files (`C12H26.xyz`). The box has the
lateral cell of the slab; its height sets the film thickness. packmol
uses its default seed, so the same command gives the same box; `--seed -1`
for a random one (independent samples). The real
thickness after the load is applied depends on the potential: check the
density in the middle of the film after the equilibration.

## 2. Stack and define the groups

```bash
python ../tools/interface.py Fe110.extxyz --middle water.extxyz --load 1.0 \
       --specorder Fe O H -o Fe_water_Fe.data
```

```
1457 atoms (432 bottom, 593 middle, 432 top), cell [19.8562 21.0607 48.3663] A -> Fe_water_Fe.data
group layer_bot   id <> 1 72
group thermo_bot  id <> 73 288
group surf_bot    id <> 289 432
group molecules   id <> 433 1025
group surf_top    id <> 1026 1169
group thermo_top  id <> 1170 1385
group layer_top   id <> 1386 1457
group sub_bot     union layer_bot thermo_bot surf_bot
group sub_top     union layer_top thermo_top surf_top
group intercalate union surf_bot molecules surf_top
group mobile      union thermo_bot intercalate thermo_top
variable    load_force equal -0.03625158341   # eV/A per atom: 1.0 GPa on 418.19 A^2, 72 atoms
-> Fe_water_Fe.groups.lmp
```

Layout, from the bottom (thicknesses from the outer surface of each slab):

| Group | Where | Default | In LAMMPS |
|---|---|---|---|
| `layer_bot` | lowest atoms of the bottom slab | `--fixed 1.5` Å | fixed |
| `thermo_bot` | above it | `--thermo 6.5` Å | thermostat |
| `surf_bot` | rest of the bottom slab | | NVE |
| `molecules` | the middle part | | NVE |
| `surf_top` | top slab near the interface | | NVE |
| `thermo_top` | | `--top-thermo` (= `--thermo`) | thermostat |
| `slide_top` | optional, `--top-slide` Å | 0 | moved (sliding) |
| `layer_top` | outermost atoms of the top slab | `--top-fixed` (= `--fixed`) | load (and sliding) |

Parts are stacked `--gap 2.5` Å apart, with `--vacuum 10` Å above the top
slab. In the LAMMPS input:

```
read_data  Fe_water_Fe.data
include    Fe_water_Fe.groups.lmp
fix        load layer_top addforce 0.0 0.0 v_load_force
```

The load: F = P·A/N, with A the xy area of the cell (418.19 Å²) and N the
atoms of `layer_top` (72): 1 GPa → 0.03625 eV/Å per atom.
`tools/load_force.py` converts in both directions, also to check an
existing input:

```bash
python ../tools/load_force.py --force -0.03625158341 --cell Fe_water_Fe.data --natoms 72
```

Checked by reading the files in LAMMPS (`read_data`, `include`,
`pair_style zero`): the group sizes are those above, `layer_bot` ends at
z = 0, `layer_top` starts at 38.37 Å, above `thermo_top` (36.34 Å).

## Variants

- **Top slab.** By default a copy of the bottom one. `--top other.extxyz`
  for a different slab with the same lateral cell (`--max-strain 0.005`):
  e.g. the film and substrate written by `zsl_interface.py` (02).
  `--flip-top` turns the top upside down, so that the same termination
  faces the interface when top and bottom terminations differ.
- **Film instead of molecules.** `--middle` takes any structure with the
  same lateral cell, e.g. an amorphous film from a melt-quench (03).
- **Type order.** With a MLIP give the full element list of the model in
  `--specorder` (e.g. `Fe C O H P Zn S`), so the same `pair_coeff` works for
  every system of the project.
- **Groups by hand.** If you build the system another way, the group
  commands are only id ranges: sort each part along z and count.
