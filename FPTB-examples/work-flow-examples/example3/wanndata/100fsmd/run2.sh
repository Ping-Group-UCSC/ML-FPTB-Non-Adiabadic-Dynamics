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
PWDIR="$HOME/q-e-qe-6.6/bin"

MPICMD="mpirun -n 1"
$MPICMD $PWDIR/wannier90.x -pp test

MPICMD="mpirun -n $SLURM_NTASKS"

$MPICMD $PWDIR/pw2wannier90.x -in pw2wann.in > pw2wann.out
