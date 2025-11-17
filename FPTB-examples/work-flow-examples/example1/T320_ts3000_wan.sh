#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH -C haswell
#SBATCH -q shared
#SBATCH -J T320_ts3000_wan.sh
#SBATCH -t 32:00:00

#OpenMP settings:
#export OMP_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

#run the application
module load espresso
module load wannier90
srun -n 1 wannier90.x -pp T320_ts3000_map.win > T320_ts3000_wan_pre.log
srun -n 4 pw2wannier90.x < T320_ts3000_map.pw2wan > T320_ts3000_pw2wan.log
srun -n 1 wannier90.x T320_ts3000_map.win > T320_ts3000_wan.log