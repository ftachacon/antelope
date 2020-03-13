#!/bin/bash

i=$1
f=$2
dir=$3

while  [ $i -lt $f ]
do

	time -p python postprocessingAForVis0.py ${dir} $i	
	i=$[$i+1]

done

echo "all files are foundxs"
