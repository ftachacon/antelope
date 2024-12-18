#!/usr/bin/env python
import numpy as np
import os
import sys


Nsnap = 16    #32
Ntop  = 320  #160

temp = np.loadtxt("nc_it000000_iky000000.dat");

#nc_it005120_iky000100.dat

for n in range(Ntop+1):
    i = n*Nsnap
    sname="nc_it" + str("%.6d"%i) + "_iky*"    
    dname="nc_it" + str("%.6d"%i) + ".dat"
    instruction="cat " + sname + " > " + dname
    os.system( instruction )
    instruction2= "rm " + sname
    os.system( instruction2 )
    print("\nIndex = ",i,";                   Source-File = ", sname)

print ("\n***************************")
print("\nEND OF CONCATENATIONS\n")
print ("\n***************************")
