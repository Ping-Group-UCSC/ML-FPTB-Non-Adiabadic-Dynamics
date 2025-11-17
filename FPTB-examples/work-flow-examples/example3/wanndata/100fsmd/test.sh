#!/bin/bash
#SBATCH --job-name=qe
#SBATCH --output=qe.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=cpuq
#SBATCH --account=cpuq
# SBATCH --dependency=afterany:

source ~/.bashrc

./rya.bash 
