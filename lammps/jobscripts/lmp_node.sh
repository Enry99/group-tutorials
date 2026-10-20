#!/bin/bash
# LAMMPS (Kokkos GPU build with ML-IAP) on Leonardo Booster: N full nodes,
# 4 MPI ranks per node, one per GPU, 8 cores each.
# For large systems (10^5 atoms and more). Not tested on the examples of
# this folder: check the speed against lmp_1gpu.sh on your system.
#
# Usage:   sbatch lmp_node.sh input.in        -> writes input.out and input.log
#
# The --qos=boost_qos_dbg line selects the debug queue (max 30 min, 2 jobs at
# a time). Delete it for production runs.
#SBATCH --job-name=lmp
#SBATCH --account=CHANGE_ME          # project budget; list yours with: saldo -b
#SBATCH --partition=boost_usr_prod   # Booster partition (4x A100 64 GB per node)
#SBATCH --qos=boost_qos_dbg          # debug queue (see above)
#SBATCH --time=00:30:00              # walltime (max 24 h without the debug qos)
#SBATCH --nodes=1                    # number of nodes: 4 GPUs each
#SBATCH --ntasks-per-node=4          # 4 MPI ranks per node ...
#SBATCH --gres=gpu:4                 # ... one per GPU
#SBATCH --cpus-per-task=8            # 8 cores per rank: the whole node
#SBATCH --mem=0                      # all the memory of the node
#SBATCH --output=%x-%j.out           # Slurm log (%x = job name, %j = job id), not the LAMMPS output
#SBATCH --error=%x-%j.err

INP=${1:?"usage: sbatch lmp_node.sh input.in"}
OUT=${INP%.in}.out
LOG=${INP%.in}.log

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load cuda/12.6 openmpi/4.1.6--gcc--12.2.0
source /leonardo/pub/userexternal/epedrett/.venv/macecuda126/bin/activate
LMP=/leonardo/pub/userexternal/epedrett/software/lammps-22Jul2026/build/lmp

# -k on g 4: each rank takes one of the 4 GPUs of its node
mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core \
    $LMP -k on g 4 -sf kk -pk kokkos newton on neigh half -in "$INP" -log "$LOG" >> "$OUT"
