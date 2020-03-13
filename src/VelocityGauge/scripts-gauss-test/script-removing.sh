#!/bin/bash

vec=(EP__-01__Phi0__--000
EP__-01__Phi0__-0.00
EP__-01__Phi0__0.06
EP__+01__Phi0__0.06
EP__-01__Phi0__0.060
EP__+01__Phi0__0.060
EP__-01__Phi0__0.2
EP__+01__Phi0__0.2
EP__-01__Phi0__0.20
EP__+01__Phi0__0.20
EP__-01__Phi0__0.40
EP__+01__Phi0__0.40
EP__-01__Phi0__0.44
EP__+01__Phi0__0.44
EP__-01__Phi0__0.5026
EP__+01__Phi0__0.5026
EP__-01__Phi0__0.57
EP__+01__Phi0__0.57
EP__-01__Phi0__0.570
EP__+01__Phi0__0.570
EP__-01__Phi0__0.75
EP__+01__Phi0__0.75
EP__-01__Phi0__1.00
EP__+01__Phi0__1.00
EP__-01__Phi0__1.000
EP__+01__Phi0__1.000
EP__-01__Phi0__1.16
EP__+01__Phi0__1.16
EP__-01__Phi0__1.160
EP__+01__Phi0__1.160
EP__-01__Phi0__.750
EP__+01__Phi0__.750
)



down=$2
Nvec=$2

for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

     echo "Dir = ${vec[iparam]}"
     j=$(( iparam +1 ))
     k=$(( iparam +2 ))
     l=$(( iparam +3 ))

     echo "i = ${iparam}"
     echo "j = $j "
     echo "k = $k "
     echo "l = $l "

     tar -cvzf "$1.tar.gz" "${vec[iparam]}"   "${vec[j]}" "${vec[k]}"   "${vec[l]}"
     
     #cd "./${vec[iparam]}"
#     rm -r c* Cher* id 
 #    rm  mgrid.dat dipoles.dat gvelocit* termial-out* exec_hh* *.sh edis*

  #   cd ../
     
     #rm HHG*dipoles.dat
     #cd ../

     echo ""
     

done
