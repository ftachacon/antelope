#!/bin/sh

#Number of MPI tasks
#SBATCH -n 23

#SBATCH -J "Wannier"

time -p mpirun -np $SLURM_NTASKS ./exec_hhg_mpi 1> stdout.log 2> stderr.log
