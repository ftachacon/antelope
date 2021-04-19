#!/bin/bash

min=$1
max=$2

for (( i=${min}; i<=${max}; i++ ))
do
    time -p python postprocessingAForVis1.py ../ $i 
    
    echo "Welcome $i times"
done
