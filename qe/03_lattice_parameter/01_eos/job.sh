#!/bin/bash
# Run every *.pwi in the current folder, one after the other, on 1 GPU.
# Used for convergence tests and E(V) scans: many short runs in one job.
#
# A run is skipped if its .pwo contains "JOB DONE." and no "convergence NOT
# achieved", so you can resubmit the same job after a timeout and it continues
# where it stopped. Crashed, killed or unconverged runs are run again.
# Careful: pw.x prints "JOB DONE." also when the SCF did not converge.
# A dry run (nstep = 0) also ends with "JOB DONE.": delete its .pwo before the real run.
#
# Usage:   sbatch job.sh            (from the folder with the .pwi files)
#SBATCH --job-name=eos-Al
#SBATCH --account=CHANGE_ME
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg          # delete for scans longer than 30 min
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load profile/chem-phys
module load quantum-espresso/7.3.1--openmpi--4.1.4--nvhpc--23.1-openblas

for INP in *.pwi; do
    OUT=${INP%.pwi}.pwo
    if [ -f "$OUT" ] && grep -q "JOB DONE." "$OUT" \
                     && ! grep -q "convergence NOT achieved" "$OUT"; then
        echo "skip $INP (done)"
        continue
    fi
    echo "run  $INP"
    mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core pw.x -nk 1 -in "$INP" > "$OUT"
done
