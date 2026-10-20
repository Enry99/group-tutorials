#!/bin/bash
# pw.x (QE 7.3.1, GPU build) on one full Leonardo Booster node: 4 GPUs x 8 cores.
# Use it for larger cells (e.g. 05: 30 atoms, many k-points).
#
# Usage:   sbatch pw_node.sh input.pwi [NPOOL]
#   NPOOL = number of k-point pools (-nk). Default 4 = one pool per GPU.
#           Use 1 for Gamma-only runs (K_POINTS gamma; 2 with nspin = 2) and when
#           there are fewer irreducible k-points than GPUs.
#SBATCH --job-name=pw
#SBATCH --account=CHANGE_ME          # project budget; list yours with: saldo -b
#SBATCH --partition=boost_usr_prod
#SBATCH --time=04:00:00              # max 24:00:00 on boost_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4          # 4 MPI ranks ...
#SBATCH --gres=gpu:4                 # ... one per GPU
#SBATCH --cpus-per-task=8            # 4 x 8 = the 32 cores of the node
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

INP=${1:?"usage: sbatch pw_node.sh input.pwi [npool]"}
NPOOL=${2:-4}
OUT=${INP%.pwi}.pwo

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load profile/chem-phys
module load quantum-espresso/7.3.1--openmpi--4.1.4--nvhpc--23.1-openblas

mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core pw.x -nk "$NPOOL" -in "$INP" > "$OUT"
