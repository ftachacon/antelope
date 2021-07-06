#!/bin/bash

module load python 

vec=(EP__-01__Phi0__0.00
EP__+01__Phi0__0.00
EP__-01__Phi0__-0.06
EP__-01__Phi0__0.06
EP__+01__Phi0__-0.06
EP__+01__Phi0__0.06
EP__-01__Phi0__-0.1
EP__-01__Phi0__0.1
EP__+01__Phi0__-0.1
EP__+01__Phi0__0.1
EP__-01__Phi0__-0.2
EP__-01__Phi0__0.2
EP__+01__Phi0__-0.2
EP__+01__Phi0__0.2
EP__-01__Phi0__0.3
EP__+01__Phi0__0.3
EP__-01__Phi0__-0.4
EP__-01__Phi0__0.4
EP__+01__Phi0__-0.4
EP__+01__Phi0__0.4
EP__-01__Phi0__0.5
EP__+01__Phi0__0.5
EP__-01__Phi0__0.6
EP__+01__Phi0__0.6
EP__-01__Phi0__0.7
EP__+01__Phi0__0.7
EP__-01__Phi0__-0.8
EP__-01__Phi0__0.8
EP__+01__Phi0__-0.8
EP__+01__Phi0__0.8
EP__-01__Phi0__0.9
EP__+01__Phi0__0.9
EP__-01__Phi0__1.00
EP__+01__Phi0__1.00
EP__-01__Phi0__1.1
EP__+01__Phi0__1.1
EP__-01__Phi0__-1.16
EP__-01__Phi0__1.16
EP__+01__Phi0__-1.16
EP__+01__Phi0__1.16
EP__-01__Phi0__1.2
EP__+01__Phi0__1.2
EP__-01__Phi0__1.3
EP__+01__Phi0__1.3
EP__-01__Phi0__1.4
EP__+01__Phi0__1.4
EP__-01__Phi0__1.5
EP__+01__Phi0__1.5
EP__-01__Phi0__-1.55
EP__-01__Phi0__1.55
EP__+01__Phi0__-1.55
EP__+01__Phi0__1.55
)



####################################
###


fname=$2
Nvec=$1 
down=$1



####################################
###
for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

    j=$(( iparam + 1 ))
    k=$(( iparam + 2 ))
    l=$(( iparam +3 ))
    echo ""
    echo "dir = ${vec[iparam]}"  
    
    tar -cvzf "${fname}.tar.gz"  "${vec[iparam]}"  "${vec[j]}"  "${vec[k]}" "${vec[l]}"
   
    echo ""
    echo ""
    echo "l = $l"
    echo "dir = ${vec[iparam]}"
    echo ""
    echo ""
 #Chern
done


echo ""
echo "\nEnd of the removing unnecesary files at ${vec[0]}\n"
echo ""
