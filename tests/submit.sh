#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 1:00:00

mpirun -np 4 a.out
