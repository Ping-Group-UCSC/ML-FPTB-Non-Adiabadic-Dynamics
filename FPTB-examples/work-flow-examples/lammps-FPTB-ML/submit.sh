#!/bin/bash
#SBATCH -A *****
#SBATCH -q regular
#SBATCH -C gpu
#SBATCH -J 20wave
#SBATCH -t 6:00:00
#SBATCH -N 8
#SBATCH -G 32
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1


#module load cpu
#module load spack
#spack load slepc/e2rxo7w
#spack load petsc/iqyxcgy

module load cudnn

module list

echo $LD_LIBRARY_PATH

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=2

export SLURM_CPU_BIND="cores"

srun -n 32 -c 32 ./ehrenfest lammps.in  > log 
