#!/bin/bash
module load python

f=$1            #Directory--name
mpitasks=$2     #Tot. No. of MPI tasks or ranks


cp postprocessing* ./$f
cd ./$f

ls -lh full*    #Verifying mpi-output files

time -p ./postprocessingBForVis1.py / 0 ${mpitasks} #Python--reduction

ls -lh int* #Verifying excitence of intra and inter currents

echo ""
echo " Reduction is finished "
echo " "

pwd 

cd ../
