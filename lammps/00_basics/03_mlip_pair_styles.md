# MLIP pair styles

The templates use MACE models through ML-IAP, the interface of the group
build. This page gives the two lines that change for each kind of model,
taken from the group inputs. In all of them `pair_coeff * *` is followed by
the element of each atom type, in the order of the data file.

## MACE through ML-IAP (templates)

```
pair_style  mliap unified MODEL.model-mliap_lammps.pt 0
pair_coeff  * * Fe O H
```

- The `.pt` file is made from the trained model, on a GPU:
  `python mace/cli/create_lammps_model.py MODEL.model --format=mliap`.
  Multi-head models must be reduced to the head you use first: ML-IAP
  takes one head.
- `0` at the end: no ghost-neighbour list for the model (as in all the
  group inputs).
- Runs with Kokkos: `-k on g N -sf kk -pk kokkos newton on neigh half`
  (`01_running_lammps.md`), one MPI rank per GPU.
- The model must contain every element of `pair_coeff`. A model trained
  on Fe, C, O, H, P, Zn, S can run a system with only Fe, O, H: list the
  types that are in the data file.

## Other MLIPs in the group inputs

These need a LAMMPS build with the corresponding package or plugin;
`lmp -h` lists the pair styles of a build. The group build is the one for
MACE through ML-IAP.

| Model | pair_style | pair_coeff |
|---|---|---|
| MACE, older interface | `mace no_domain_decomposition` (remove the keyword with more than one GPU) | `* * MODEL-lammps.pt Si C O H N` |
| MACE, symmetrix | `symmetrix/mace` | `* * MODEL.json Si C O H N` |
| DeePMD | `deepmd graph.pb` (or `compressed.pb`) | `* * Si C O H N` |
| GRACE | `grace pad_verbose` | `* * /path/to/final_model/ Si C O H N` |
| GRACE-FS | `grace/fs` | `* * FS_model.yaml Si C O H N` |
| SevenNet (parallel) | `e3gnn/parallel` | `* * 3 /path/to/deployed_parallel Si C O H N` (3 = number of message-passing layers of the model) |
| CHGNet (chgnet-lmp) | `chgnet/gpu /path/to/potentials/CHGNET` | `* * 0.3.0 C Si` (model version) |
| GAP (QUIP) | `quip` | `* * GAP.xml "Potential xml_label=..." 14 7` (atomic numbers) |

## Checks with a new model

1. **Energy conservation**: `02_md_ensembles/nve.in` at the highest
   temperature you will use. The total energy must not drift.
2. **Equilibrium structure**: lattice parameter (`01_minimize/minimize.in`
   with box/relax) and density at T (`npt.in`) against DFT or experiment.
3. **Speed**: `Performance:` line of the log, steps per second per atom.
   It sets the size and length of what you can run.

Validation against DFT (errors, parity plots) belongs to the MLIP part.
