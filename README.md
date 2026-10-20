# Group tutorials: simulations on Leonardo

Examples and templates for the codes we use: commented inputs,
Leonardo jobscripts, what to check in the output, and the errors you will
meet. Use them to learn a calculation the first time, then as a reference
when you need a template or a quick explanation.

Start from `leonardo/` (the cluster), then the folder of the code you need:

- `leonardo/`: nodes, storage, Slurm, Python environment, data transfer;
- `qe/`: Quantum ESPRESSO (pw.x, neb.x), with its own `README.md`;
- `structures/`: structure generation (slabs, interfaces, amorphous
  starting structures, sliding interfaces, tips), for any code;
- to come: LAMMPS, MLIP training and active learning.

## Setup (once)

Clone the repository in your area and work in a copy of the code folder,
so that the clone stays clean and can be updated with `git pull`:

```bash
git clone <group-repo-url> group-tutorials
cp -r group-tutorials/qe qe_runs
```

Then follow the setup in the code's `README.md` (for QE: your project in
the jobscripts, Python aliases).

## Found a problem?

1. Search the code's `troubleshooting.md` for the exact error text.
2. Not there? Open an issue in the group repository.
3. Once solved, add it to `troubleshooting.md`: message, cause, fix.
