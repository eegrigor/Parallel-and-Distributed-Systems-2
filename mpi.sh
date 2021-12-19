#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=5:00

module load gcc openmpi

mpicc retroactive.c -lm -o mpi

srun ./mpi