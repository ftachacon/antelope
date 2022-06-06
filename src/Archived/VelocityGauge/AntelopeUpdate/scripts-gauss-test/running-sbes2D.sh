#!/bin/bash
#SBATCH -t 3:59:00              #Walltime
##SBATCH -N 3                   #Number of nodes (commented)
#SBATCH -n 119                  #Number of MPI tasks
##SBATCH --ntasks-per-node=16   #Number of MPI tasks or ranks per node (commented)
#SBATCH -o ChernI-HM.log
#SBATCH -e ChernI-HM.err
#SBATCH --qos=interactive             #standard    
#SBATCH --job-name=ChernIns
#SBATCH --mail-user=ftachacon@gmail.com
#SBATCH --mail-typ=ALL

source /etc/bashrc
source ~/.bash_profile

#What every you did to compile code, do here
module load mkl openmpi intel
module list

cd $SLURM_SUBMIT_DIR

#DEBUG
ls -l $EX_DIR/exec_hhg_mpi
#DEBUG -- while we are at it, verify:
echo "NTASKS=" $SLURM_NTASKS

pwd
date
mpirun -n $SLURM_NTASKS exec_hhg_mpi $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} > termial-outpput.dat
date
