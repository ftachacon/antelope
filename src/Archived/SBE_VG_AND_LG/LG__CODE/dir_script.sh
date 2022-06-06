#!/bin/bash

#e=$1
#p=$2
#mt2=$3
#dir="EP__${e}__Phi__${p}__Mt2__${mt2}"

dir=$1
mkdir ${dir}

cp exec* ${dir}/
cp mpi_gcc-4.9.2.sh ${dir}
cp input.dat ${dir}

