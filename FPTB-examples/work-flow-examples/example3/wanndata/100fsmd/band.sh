#!/bin/bash
#SBATCH --job-name=qe
#SBATCH --output=qe.%j.out
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=24:00:00
#SBATCH --partition=cpuq
#SBATCH --account=cpuq
# SBATCH --dependency=afterany:

module load intel/impi
MPICMD="mpirun -n $SLURM_NTASKS "
PWDIR="$HOME/q-e-qe-6.6/bin"

#$MPICMD $PWDIR/pw.x -nk 4 -nd 144  -inp scf.in > scf.log
$MPICMD $PWDIR/pw.x -nk 4 -nd 144  -inp bands.in > bands.log
