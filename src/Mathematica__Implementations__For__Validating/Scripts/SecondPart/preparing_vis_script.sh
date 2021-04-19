#!/bin/bash

i=$1
f=$2


while  [ $i -lt $f ]
do

	time -p python postprocessingAForVis0.py ../ $i	
	i=$[$i+1]

done

echo "all files are foundxs"
