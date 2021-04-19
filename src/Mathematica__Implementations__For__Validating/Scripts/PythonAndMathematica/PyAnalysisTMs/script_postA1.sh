#!/bin/bash

min=$1
max=$2
dir=$3

cp ../outlaserdata.dat .

for (( i=${min}; i<=${max}; i++ ))
do
    time -p python postprocessingAForVis1.py ${dir} $i 
    
    echo "Welcome $i times"
done
