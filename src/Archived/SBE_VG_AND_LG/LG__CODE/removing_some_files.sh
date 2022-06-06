#!/bin/bash
dir=$1

cd ./${dir}

file="interband_dipole_full_evol.dat"

if [ -f "${file}" ]
then

    echo '${file} found. then, full_* will be removed. '
    rm -r -f full_integrated_currents_rank_0*    
    rm -r -f connection.dat curvature.dat gvelocities.dat dipoles.dat edispersion.dat mgrid.dat

    echo "Removing of full_integ* is done"
else

    echo "${file} not found."

fi
 
 
echo ""

cd ..
 
