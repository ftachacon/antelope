#!/bin/bash

#module load python 

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

mpi=119
Nvec=$2 
down=$1
file="intraband_current_full_evol.dat"


####################################
###
for((iparam=$down; iparam<=${Nvec}; iparam++ ))
do 

    echo ""
    echo "dir = ${vec[iparam]}"  
    
    cd ./"${vec[iparam]}"
    
    #"./"reduction.sh  "${vec[iparam]}" ${mpi}
    #"./"removing_mpi_files.sh  "${vec[iparam]}"   
    
    rm -r HHG*  exec_hhg_mpi
    rm *.py 
    rm -r *Figure launchingJobs.sh~
    rm *.sh  id 
    rm Chern*
    cd ..
    #rm 
     echo ""
     echo "dir = ${vec[iparam]}"    
 #Chern
done


echo ""
echo "\nEnd of the removing unnecesary files at ${vec[0]}\n"
echo ""

