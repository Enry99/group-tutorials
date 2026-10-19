#!/bin/bash
# neb.x (QE 7.3.1, GPU build) on one Leonardo Booster node: 4 image groups x 1 GPU.
#
# The first and last images are fixed: after the first step only the
# num_of_images - 2 intermediate images are computed. Make that number equal
# to the number of image groups (-ni), so every GPU always has an image:
#   num_of_images = 6  -> 4 moving images -> 1 node,  -ni 4  (this script)
#   num_of_images = 10 -> 8 moving images -> 2 nodes, -ni 8  (--nodes=2)
# The number of MPI ranks must be a multiple of -ni.
#
# Usage:   sbatch job.sh neb.pwi
#SBATCH --job-name=neb-H
#SBATCH --account=CHANGE_ME
#SBATCH --partition=boost_usr_prod
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

INP=${1:?"usage: sbatch job.sh neb.pwi"}
OUT=${INP%.pwi}.pwo

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load profile/chem-phys
module load quantum-espresso/7.3.1--openmpi--4.1.4--nvhpc--23.1-openblas

# -ni = one image group per MPI rank (1 GPU each)
mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core neb.x -ni $SLURM_NTASKS -input "$INP" > "$OUT"
