#!/bin/bash
# LAMMPS (Kokkos GPU build with ML-IAP) on Leonardo Booster: 1 GPU + 8 cores.
# Enough for MLIP runs up to ~10^4-10^5 atoms (depends on the model).
#
# Usage:   sbatch lmp_1gpu.sh input.in        -> writes input.out and input.log
#
# Lines starting with "#SBATCH" are read by Slurm before the job starts.
# The --qos=boost_qos_dbg line selects the debug queue (max 30 min, 2 jobs at
# a time): good for tests and tutorials. Delete it for production runs.
#SBATCH --job-name=lmp
#SBATCH --account=CHANGE_ME          # project budget; list yours with: saldo -b
#SBATCH --partition=boost_usr_prod   # Booster partition (4x A100 64 GB per node)
#SBATCH --qos=boost_qos_dbg          # debug queue (see above)
#SBATCH --time=00:30:00              # walltime (max 24 h without the debug qos)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1          # 1 MPI rank ...
#SBATCH --gres=gpu:1                 # ... driving 1 GPU
#SBATCH --cpus-per-task=8            # 8 cores per rank: a quarter of a node
#SBATCH --mem=100000                 # MB, a quarter of the node memory
#SBATCH --output=%x-%j.out           # Slurm log (%x = job name, %j = job id), not the LAMMPS output
#SBATCH --error=%x-%j.err

# Input file from the command line; screen output -> .out, LAMMPS log -> .log
INP=${1:?"usage: sbatch lmp_1gpu.sh input.in"}
OUT=${INP%.in}.out
LOG=${INP%.in}.log

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Group build of LAMMPS (stable 22 Jul 2025, Kokkos CUDA + ML-IAP for MACE).
# It needs CUDA, OpenMPI and the shared Python environment (torch, mace).
module purge
module load cuda/12.6 openmpi/4.1.6--gcc--12.2.0
source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate
LMP=/leonardo/pub/userexternal/epedrett/software/lammps-22Jul2026/build/lmp

# -k on g 1: Kokkos on 1 GPU; -sf kk: Kokkos versions of the styles;
# -pk kokkos newton on neigh half: the settings ML-IAP needs
mpirun -np 1 $LMP -k on g 1 -sf kk -pk kokkos newton on neigh half -in "$INP" -log "$LOG" >> "$OUT"
