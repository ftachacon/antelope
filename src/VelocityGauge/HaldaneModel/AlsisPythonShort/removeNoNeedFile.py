#!/usr/bin/env python                                                                                                   
import numpy as np
import os
import sys

Nsnap = 32
Ntop  = 160
i=32

for n in range(Ntop+1):
    i=n*Nsnap
    sname0  ="nc_it" + str("%.6d"%i) + "_iky000101.dat"
    sname1 ="nc_it" + str("%.6d"%i) + "_iky000102.dat"
    instruction0="rm " + sname0
    instruction1="rm " + sname1
    os.system( instruction0 )
    os.system( instruction1 )
