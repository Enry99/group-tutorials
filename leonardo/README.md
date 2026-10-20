# Working on Leonardo (CINECA)

What every code needs: the nodes, where to put files, how to submit jobs,
the Python environment and how to move data. The code-specific parts
(modules, jobscripts, checks) are in each code's folder, e.g.
`qe/00_basics/01_running_qe.md`. More info on the cluster at [CINECA documentation](https://docs.hpc.cineca.it).

## Logging in

CINECA uses two-factor authentication. Once: activate it at
https://sso.hpc.cineca.it (password and OTP app) and install the
`smallstep` client (`step`; not the unrelated `step` package of the Linux
repositories), configured as in the CINECA page
[Access to the systems](https://docs.hpc.cineca.it/general/access.html).
Then, every 12 hours (the certificate expires):

```bash
step ssh login '<your-email>' --provisioner cineca-hpc   # opens the browser: password + OTP
ssh <username>@login.leonardo.cineca.it
```

## The machine

- **Booster** nodes: 32 CPU cores and 4 NVIDIA A100 GPUs (64 GB) each.
  Our codes (QE, LAMMPS, MLIP training) run on the GPUs.
- **Login nodes** are shared, with a 10-minute CPU limit per process: edit,
  submit, run short Python scripts. Never run simulations there.
- **Budget**: `saldo -b` lists your projects and the hours left; the project
  name goes in `#SBATCH --account=` of the jobscripts.

## Where to put files

| Area | Use | Note |
|---|---|---|
| `$HOME` | scripts, small files, small Python environments | 50 GB, backed up |
| `$WORK`, `$FAST` | runs, results | per project, no backup, deleted 6 months after the project ends |
| `$SCRATCH` | large temporary files | files not accessed for 40 days are deleted |

## Slurm

A job is a bash script (`job.sh`) whose `#SBATCH` lines ask for resources;
Slurm runs it on a compute node when they are free.

```bash
sbatch job.sh <arguments>   # submit; prints the job id
squeue --me                 # your jobs: PD = pending, R = running
scancel <jobid>             # kill a job
saldo -b                    # remaining budget
```

The single-GPU and loop jobscripts use the debug queue (`#SBATCH --qos=boost_qos_dbg`:
at most 30 min and 2 jobs at a time), good for tests and tutorials. For
production runs delete that line and set `--time` (at most 24 h on
`boost_usr_prod`).

## Python environment

The scripts of the tutorials need Python with a few packages (ase, numpy,
scipy, pandas, spglib, matplotlib; pymatgen for `structures/`, OVITO for
part of `lammps/analysis/`, which is in the environment
`/leonardo/pub/userexternal/epedrett/.venv/mace`). First try the shared
environment:

```bash
source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate
```

To load it at every login, add that line to your `~/.bashrc`:

```bash
echo "source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate" >> ~/.bashrc
```

To create your own environment instead (`YOUR_ENV_NAME` is a name of your
choice):

```bash
module av python                 # pick Python >= 3.10
module load python
mkdir -p $HOME/venvs
python -m venv $HOME/venvs/YOUR_ENV_NAME
source $HOME/venvs/YOUR_ENV_NAME/bin/activate
pip install ase numpy scipy pandas spglib matplotlib pymatgen   # ovito if needed
echo "source $HOME/venvs/YOUR_ENV_NAME/bin/activate" >> ~/.bashrc
```

## Moving data

For large transfers use the data mover `data.leonardo.cineca.it` (no
shell, no CPU limit) with absolute paths. From your computer:

```bash
rsync -rav <user>@data.leonardo.cineca.it:/absolute/path/to/remote/folder/ ./local_folder/
```

or, on Windows, scp:

```bash
scp -r <user>@data.leonardo.cineca.it:/absolute/path/to/remote/folder/ ./local_folder/
```
