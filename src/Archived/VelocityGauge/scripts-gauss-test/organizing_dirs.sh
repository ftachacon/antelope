#!/bin/bash 

vec=(EP__-01__Phi0__0.00
EP__+01__Phi0__0.00
EP__-01__Phi0__0.06283185
EP__+01__Phi0__0.06283185
EP__-01__Phi0__0.12566371
EP__+01__Phi0__0.12566371
EP__-01__Phi0__0.18849556
EP__+01__Phi0__0.18849556
EP__-01__Phi0__0.25132741
EP__+01__Phi0__0.25132741
EP__-01__Phi0__0.31415927
EP__+01__Phi0__0.31415927
EP__-01__Phi0__0.37699112
EP__+01__Phi0__0.37699112
EP__-01__Phi0__0.43982297
EP__+01__Phi0__0.43982297
EP__-01__Phi0__0.50265482
EP__+01__Phi0__0.50265482
EP__-01__Phi0__0.56548668
EP__+01__Phi0__0.56548668
EP__-01__Phi0__0.62831853
EP__+01__Phi0__0.62831853
EP__-01__Phi0__0.69115038
EP__+01__Phi0__0.69115038
EP__-01__Phi0__0.75398224
EP__+01__Phi0__0.75398224
EP__-01__Phi0__0.81681409
EP__+01__Phi0__0.81681409
EP__-01__Phi0__0.87964594
EP__+01__Phi0__0.87964594
EP__-01__Phi0__0.9424778
EP__+01__Phi0__0.9424778
EP__-01__Phi0__1.00530965
EP__+01__Phi0__1.00530965
EP__-01__Phi0__1.0681415
EP__+01__Phi0__1.0681415
EP__-01__Phi0__1.13097336
EP__+01__Phi0__1.13097336
EP__-01__Phi0__1.19380521
EP__+01__Phi0__1.19380521
EP__-01__Phi0__1.25663706
EP__+01__Phi0__1.25663706
EP__-01__Phi0__1.31946891
EP__+01__Phi0__1.31946891
EP__-01__Phi0__1.38230077
EP__+01__Phi0__1.38230077
EP__-01__Phi0__1.44513262
EP__+01__Phi0__1.44513262
EP__-01__Phi0__1.50796447
EP__+01__Phi0__1.50796447
)

down=3
Nvec=49



####################################
###
for((iparam=${down}; iparam<=${Nvec}; iparam++ ))
do 

     echo ""
     echo "dir = ${vec[iparam]}"  
    
     mv "${vec[iparam]}" ./Phi_Half_Region

    #"./"reduction.sh  "${vec[iparam]}" 
    #"./"removing_mpi_files.sh  "${vec[iparam]}"   

     echo ""
     echo "dir = ${vec[iparam]}"    

done
