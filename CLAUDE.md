# Working on this repository

Reference for the thesis students of the group (E. Pedretti): templates,
examples and troubleshooting for Quantum ESPRESSO, LAMMPS with MLIPs and
structure generation on CINECA Leonardo. Students copy a template and
adapt it; it is not a step-by-step course. Accessibility comes before
completeness.

## Layout

- `leonardo/`: login, storage, Slurm, Python environments.
- `qe/`: numbered examples, `jobscripts/`, `tools/`, `troubleshooting.md`.
  Checked by the owner: change it only for real errors and list every
  change in the reply.
- `structures/`: `tools/` (slab, interface, density, packmol, tip, ...)
  and numbered examples; builds the data and groups files of `lammps/`.
- `lammps/`: `00_basics/`, templates `01`-`09`, `jobscripts/`,
  `troubleshooting.md`, `analysis/` (tools, `tests/synthetic.py`).
- `TODO.md`: open points; keep it current.

## Writing

- Repository in English; chat with the owner in Italian.
- Plain, direct prose, short sentences, no mannered phrasing, not verbose.
  Say what to do and what to check; cut process details a student does not
  need (they go in `reference/`).
- One explanation per topic, linked from elsewhere (troubleshooting holds
  the error messages).
- Every number comes from a run in this repository or from a cited
  source. References: look up each DOI (Crossref/DataCite) and check title,
  first author, year; never write a DOI, volume or page from memory.
- No model names or tool attributions inside repository files.

## Conventions of the code

- Owner's scripts come first: when a tool derives from one, keep its
  choices and list the differences in the README ("Differences from the
  original scripts"); ask before changing a choice that is not a bug.
- Discuss before changes to the workflow (new stages, defaults such as the
  timestep, file layout).
- LAMMPS templates: settings block at the top, times in ps (steps from
  `$(round(v_t/dt))`), `units metal`, `CHANGE_ME` for account and model,
  `variable elements` in type order, `dump custom ... id element xu yu zu`
  with `sort id`, one `write_restart` per stage; an input that starts a new
  stage from a restart begins with `reset_timestep 0` (the analysis tools
  split stages on it). Slabs and interfaces: `processors * * 1` and the
  x, y balance.
- Jobscripts: each folder has `job.sh`, a copy of `jobscripts/` that
  differs only in the usage line or job name; edit them together.
- Python tools: docstring with Usage / What it does / Needs, argparse,
  `-h`; outputs `<out>.csv` and `<out>.png` (matplotlib Agg, ticks in);
  shared helpers in `common.py`, `lmp_io.py`, `ovito_io.py`.

## Checks before committing

- Run what you changed. LAMMPS templates: a local check with a classical
  potential in place of the MLIP (CPU build), reduced settings stated in
  the README; logs and figures in `reference/` as `*_cpu-serial.*`.
- `structures/tools`: the data and groups files of `lammps/05`-`08` must
  rebuild byte for byte from the commands in their READMEs.
- `lammps/analysis`: `python tests/synthetic.py` (all checks pass);
  regenerate `reference/commands_output.txt` when the output changes.
- `pyflakes` on the Python files; paths named in the READMEs exist.

## Git

- Work on a branch, push it, open a pull request only when asked.
- Commit messages: a summary line, then a short list of what changed.

## Leonardo paths (for the documentation)

- Shared Python: `/leonardo/pub/userexternal/epedrett/.venv/macecuda126`;
  OVITO: `/leonardo/pub/userexternal/epedrett/.venv/mace`.
- LAMMPS group build (30 Mar 2026, 16f8568; the directory name is not the
  version): `/leonardo/pub/userexternal/epedrett/software/lammps-22Jul2026/build/lmp`.
- packmol: `/leonardo/pub/userexternal/epedrett/software/packmol-21.2.1`;
  pseudopotentials: `/leonardo/pub/userexternal/epedrett/pseudo`.
