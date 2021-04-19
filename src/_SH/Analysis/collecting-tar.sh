#!/bin/bash

down=$1
Nvec=$2

####################################
###
for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

    j=$(( iparam*4 ))
    
    echo ""
    echo "my-j = $j"
    echo "my-i = ${iparam}\n"
    
    ./tar-compresing.sh ${j} "00${iparam}" 
   
    echo ""  
    echo "my-j = $j"

    echo ""

done
