#!/bin/bash
# Timestep test on Leonardo Booster: one Slurm array task per timestep, each
# on 1 GPU (as lmp_1gpu.sh), all from the same starting state.
#
# Usage:   sbatch job_scan.sh nve.in                 -> nve_dt0.25.log, nve_dt0.5.log, ...
#          sbatch job_scan.sh layout.in 800          -> layout_dt*.log, thermostats at 800 K
#          sbatch --array=3 job_scan.sh nve.in 300 2.014   -> 1 fs with deuterium mass
# nve.in starts from equil.restart: submit it after equil.in,
#          jid=$(sbatch --parsable job.sh equil.in)
#          sbatch --dependency=afterok:$jid job_scan.sh nve.in
# Then:    python check_dt.py "nve_dt*.log"
#SBATCH --job-name=dtscan
#SBATCH --account=CHANGE_ME          # project budget; list yours with: saldo -b
#SBATCH --partition=boost_usr_prod   # Booster partition (4x A100 64 GB per node)
#SBATCH --time=04:00:00              # the smallest dt takes longest: t_run / dt steps
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100000
#SBATCH --array=0-4                  # one task per entry of DTS
#SBATCH --output=%x-%A_%a.out        # %A = array job id, %a = task index
#SBATCH --error=%x-%A_%a.err

DTS=(0.25 0.5 0.75 1.0 1.25)         # fs; the smallest is the reference of check_dt.py
INP=${1:?"usage: sbatch job_scan.sh nve.in|layout.in [T] [mass of H]"}
T=${2:-300}                          # K, used by layout.in (nve.in keeps the T of equil.in)
MH=${3:-0}                           # 0: mass of H as in the data file
DT=${DTS[$SLURM_ARRAY_TASK_ID]}
TAG=${INP%.in}_dt${DT}
[ "$MH" != 0 ] && TAG=${TAG}_mH${MH}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Group build of LAMMPS (30 Mar 2026, development version 16f8568; Kokkos CUDA + ML-IAP for MACE;
# the directory name is not the version).
module purge
module load cuda/12.6 openmpi/4.1.6--gcc--12.2.0
source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate
LMP=/leonardo/pub/userexternal/epedrett/software/lammps-22Jul2026/build/lmp

# --map-by socket:PE=N: the rank gets its N cores (by default OpenMPI binds
# one rank to one core, so all its threads would share it)
mpirun -np 1 --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core \
    $LMP -k on g 1 -sf kk -pk kokkos newton on neigh half -in "$INP" \
    -var dt "$DT" -var T "$T" -var mH "$MH" -var tag "$TAG" -log "$TAG.log" >> "$TAG.out"
