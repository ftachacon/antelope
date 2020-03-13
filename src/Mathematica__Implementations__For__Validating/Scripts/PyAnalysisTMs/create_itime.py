#!/usr/bin/env python
import numpy as np
import os 
import sys


N=160
snap=32

x=np.zeros(N);

for n in range(N):
    x[n] = n*snap

np.savetxt("timeindex.dat",x,"%.6d")
