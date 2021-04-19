#!/bin/bash 
module load python

dir=$1 #Directory where info and scripts will be copied


cp postprocessingBForVis1.py ./${dir}
cp exec_hhg_mpi ./${dir}
cp *.sh ./${dir}


cd ./${dir}
ls -lh 

cd ..

echo "${dir}"
echo 
