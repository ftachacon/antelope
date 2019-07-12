#!/bin/bash

N_MPI_TASKs=1 #No. of MPI task or cores for parallelization



###################################
mkdir Test #creating test directory



####################################
make #compiling code and generating executable 



####################################
#Copying exe* and other files to Test dir.
cp exec_hhg_mpi ./Test/
cp -R AlsisPythonShort ./Test/
cp AlsisPythonShort/postprocessingBForVis1.py ./Test/



###################################
#Running test example
cd ./Test

date

time -p mpirun -n ${N_MPI_TASKs} exec_hhg_mpi 0. 2.54 100 11 0.006 5 1.25 440 0.0 0.0 2

date 


###################################
#Next line performe remain
echo "--------------------"
echo ""
echo "Integration along ky direction via python script"

time -p python postprocessingBForVis1.py / 0 ${N_MPI_TASKs}

echo ""
echo "End of the integration along ky-direction"


###################################
##Analysis-first-step 
cd  ./AlsisPythonShort/



######################################
##Convering data to python analysis codee
./script_to_analysis.sh 



cd ../
echo ""
echo "" 
