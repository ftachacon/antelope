#!/bin/bash

module load python 

vec=( #Set of directories
)


####################################
###

down=$1
Nvec=$2
MPITasks=$3


file="intraband_current_full_evol.dat"


####################################
###
for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

     echo ""
     echo "dir = ${vec[iparam]}"  
    
    "./"reduction.sh "${vec[iparam]}"  ${MPITasks}
    "./"removing_mpi_files.sh  "${vec[iparam]}"   

     echo ""
done


echo ""
echo "\nEnd of reduction and removing files\n"
echo ""

