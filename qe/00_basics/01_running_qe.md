# Running Quantum ESPRESSO on Leonardo

How QE runs on the Booster nodes, which jobscript to use, what to check in
a GPU run, and the Python tools to look at QE files. Nodes, storage, Slurm
and the Python environment: `../../leonardo/README.md`.

## QE on the GPUs

pw.x runs with one MPI rank per GPU and 8 OpenMP threads per rank, i.e. a
quarter of a node per rank. The GPU build is loaded by every jobscript:

```bash
module load profile/chem-phys
module load quantum-espresso/7.3.1--openmpi--4.1.4--nvhpc--23.1-openblas
```

`module av quantum-espresso` lists the versions. Do not mix versions
within a project.

## The jobscripts

In `jobscripts/`, copied into each example as `job.sh`. Each line is
commented: read one once from top to bottom.

| Script | Resources | Used for |
|---|---|---|
| `pw_1gpu.sh` | 1 GPU, 8 cores, debug queue (30 min, max 2 jobs) | bulk, slabs and molecules of the examples |
| `pw_loop.sh` | same, runs every `*.pwi` of the folder | convergence scans, E(V) |
| `pw_node.sh` | 4 GPUs, 32 cores | larger cells (05) |
| `neb.sh` | 4 GPUs, one per NEB image | neb.x (08) |

For production runs: delete the `--qos=boost_qos_dbg` line, adjust
`--time` in the jobscript and set `max_seconds` in the `.pwi` accordingly.

## Checks on a GPU run

- Near the top of the `.pwo`: `GPU acceleration is ACTIVE.`
- k-point pools (`-nk N` in the jobscript): the k-points are split into N
  groups, one per GPU. Use N = number of GPUs when there are at least as
  many irreducible k-points. With `K_POINTS gamma` you can use `-nk 2` for
  spin-polarized calculations.
- Out of GPU memory (the job dies early, error in the Slurm `.err` file):
  see "Memory errors" in `../troubleshooting.md`.

## Python tools

The scripts in `tools/` need the Python environment of
`../../leonardo/README.md`. Aliases for the three you will use often (add
them to `~/.bashrc`; with your own environment, put its `python` in place
of the shared one):

```bash
alias qe2xyz='/leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/python /leonardo/pub/userexternal/epedrett/tools/qe/qe2xyz.py'
alias pwosum='/leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/python /leonardo/pub/userexternal/epedrett/tools/qe/pwo_summary.py'
alias orthocell='/leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/python /leonardo/pub/userexternal/epedrett/tools/qe/ortho_cell.py'
```

## Looking at structures

QE files are poorly supported by viewers: XCrySDen reads `.pwi`/`.pwo` but
is old and Linux-only, OVITO reads `.pwo` only in the Pro version, recent
ASE reads `.pwo` but not `.pwi` with `ibrav ≠ 0`. The simplest route:
convert to extended XYZ and open it in OVITO (basic), ASE gui or any
other viewer.

```bash
qe2xyz slab.pwi            # any ibrav, also NEB inputs -> slab.extxyz
qe2xyz relax.pwo           # all ionic steps; add --last for the final one only
qe2xyz neb.axsf            # all images of a NEB path
qe2xyz Fe2O3.pwi -f cif    # other formats, e.g. CIF for VESTA
```

The `.extxyz` of an input also contains the QE species (`qe_type`: colour
by it to see Fe1/Fe2), the fixed atoms (`fixed`) and the starting
magnetization (`start_mag`). On your laptop: `pip install ase numpy`, copy
`tools/qe2xyz.py`, then `python qe2xyz.py slab.pwi` and so on.
