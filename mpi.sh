#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=5:00

module load gcc openmpi

mpicc distribute.c -lm -o mpi

srun ./mpi
