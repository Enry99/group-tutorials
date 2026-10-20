# Group tutorials: simulations on Leonardo

Templates and examples for the codes we use: commented inputs, Leonardo
jobscripts, what to check in the output, and the errors you will meet. It
is a reference rather than a course: find the example closest to your
calculation, copy it, change the settings at the top, and check the
output as its README says. Errors: `troubleshooting.md` of each code.

Start from `leonardo/` (the cluster), then the folder of the code you need:

- `leonardo/`: logging in, nodes, storage, Slurm, Python environment, data
  transfer;
- `qe/`: Quantum ESPRESSO (pw.x, neb.x), with its own `README.md`;
- `structures/`: structure generation (slabs, interfaces, amorphous
  starting structures, sliding interfaces, tips), for any code;
- `lammps/`: LAMMPS templates (minimization, NVE/NVT/NPT, melt-quench,
  deformation, sliding, indentation, reactions, deposition) with MLIPs,
  and the analysis of their outputs (`lammps/analysis/`);
- to come: MLIP training and active learning.

## Setup (once)

Clone the repository in your project area (`$WORK`, see
`leonardo/README.md`; `$HOME` is small) and work in a copy of the code
folder, so that the clone stays clean and can be updated with `git pull`:

```bash
git clone https://github.com/Enry99/group-tutorials.git group-tutorials
cp -r group-tutorials/qe qe_runs
```

Then follow the setup in the code's `README.md` (for QE: your project in
the jobscripts, Python aliases).

## Found a problem?

1. Search the code's `troubleshooting.md` for the exact error text.
2. Not there? Open an issue in the group repository.
3. Once solved, add it to `troubleshooting.md`: message, cause, fix.
