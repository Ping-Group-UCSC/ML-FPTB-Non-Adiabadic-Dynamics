#!/bin/bash
#SBATCH -N 4
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J T320_ts3000_dft.sh
#SBATCH -t 8:00:00

#OpenMP settings:
#export OMP_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

#run the application
module load espresso
srun -n 64 -c 4 --cpu_bind=cores pw.x -i T320_ts3000_scf.in -npool 8 > T320_ts3000_scf.log
srun -n 64 -c 4 --cpu_bind=cores pw.x -i T320_ts3000_nscf.in -npool 8 > T320_ts3000_nscf.log