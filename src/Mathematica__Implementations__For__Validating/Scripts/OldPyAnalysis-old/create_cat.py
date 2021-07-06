#!/usr/bin/env python
import numpy as np
import os
import sys


Nsnap = 32
Ntop  = 160

#nc_it005120_iky000100.dat

for n in range(Ntop+1):
    i = n*Nsnap
    sname="nc_it" + str("%.6d"%i) + "_iky*"
    dname="nc_it" + str("%.6d"%i) + ".dat"
    instruction="cat " + sname + " > " + dname
    os.system( instruction )
    instruction2= "rm " + sname
    os.system( instruction2 )

print("\nEND OF CONCATENATIONS\n")