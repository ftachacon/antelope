#!/bin/bash

i=$1 #inital index limit for ky
f=$2 #final index limit for ky points


while  [ $i -lt $f ]
do

	time -p python postprocessingAForVis0.py ../ $i	
	i=$[$i+1]

done

echo "all files are foundxs"
