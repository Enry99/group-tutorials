#!/bin/bash
# pw.x (Quantum ESPRESSO 7.3.1, GPU build) on Leonardo Booster: 1 GPU + 8 cores.
# Use it for the small cells of the examples (bulk, Al slabs, molecules in a box).
#
# Usage:   sbatch job.sh input.pwi        -> writes input.pwo
#
# Lines starting with "#SBATCH" are read by Slurm before the job starts.
# They must come before the first command of the script.
# The --qos=boost_qos_dbg line selects the debug queue (max 30 min, 2 jobs at
# a time): good for tests and tutorials. Delete it for production runs.
#SBATCH --job-name=Al111-H-fcc2
#SBATCH --account=CHANGE_ME          # project budget; list yours with: saldo -b
#SBATCH --partition=boost_usr_prod   # Booster partition (4x A100 64 GB per node)
#SBATCH --qos=boost_qos_dbg          # debug queue (see above)
#SBATCH --time=00:30:00              # walltime. Set max_seconds in the input ~5% lower
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1          # 1 MPI rank ...
#SBATCH --gres=gpu:1                 # ... driving 1 GPU
#SBATCH --cpus-per-task=8            # 8 cores per rank: a quarter of a node
#SBATCH --output=%x-%j.out           # Slurm log (%x = job name, %j = job id), not the QE output
#SBATCH --error=%x-%j.err

# Input file from the command line; the output gets the same name with .pwo
INP=${1:?"usage: sbatch job.sh input.pwi"}
OUT=${INP%.pwi}.pwo

# OpenMP threads per MPI rank = cores per rank, pinned to physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Start from a clean environment, then load the QE GPU build
module purge
module load profile/chem-phys
module load quantum-espresso/7.3.1--openmpi--4.1.4--nvhpc--23.1-openblas

mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core pw.x -nk 1 -in "$INP" > "$OUT"
