#!/bin/bash

awk '{ print $4 , $5 }' full_integrated_currents.dat > intraband_current_full_evol.dat
awk '{ print $2 , $3 }' full_integrated_currents.dat >interband_dipole_full_evol.dat

