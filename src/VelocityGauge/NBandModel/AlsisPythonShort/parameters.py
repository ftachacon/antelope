#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
pi      = np.pi;
I       = 1j;


###############################
Nvects      = 3;
Ndim        = 2;
Nsnap       = 160; #snaptshot per opt. cycle


###############################
#Lattice constant
a0          = 1./.529;
t1          = .075;
t2          = t1/3.;
M0          = 2.54*t2;
phi0        = 1.16;


###############################
## MOMENTUM AXIS
Nx      = 401;    #201
Ny      = 101;    #101; #101; #;


###############################
## LASER
I0      = 1.26e12;
w0      = 0.014;
cep0    = 0.;
Ncycles = 8.;
dt      = 1.25;
Ntime   = 15245


