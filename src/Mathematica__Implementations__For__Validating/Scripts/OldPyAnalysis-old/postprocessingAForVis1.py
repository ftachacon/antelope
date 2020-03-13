#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:53:15 2019

@author: achacon
"""

import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator 
import numpy as np
import os 
import sys 
import shutil 
import cmath
import errno

from scipy.fftpack import fft, ifft
from scipy.integrate import odeint
from scipy.integrate import ode

import time
from scipy import special
from scipy.integrate import quadrature, quad
import cstructure as cs
#import myfuncs


################################
plt.rc('axes',lw=2.2)


## CONSTANTS
pi      = np.pi;
Nsnap   = 160; #snaptshot per opt. cycle


################################
### INPUT PARAMETS
set_of_params   = sys.argv;


iparam          = set_of_params[1]
jparam          = int(set_of_params[2])


################################
##  Basic current path
BasicPath       = os.getcwd();





################################
#Need to load files
# OccupationCBMath_ncNx_00250_kyIndex_
# CoherenceMathRealPiNx_00250_kyIndex_
FolderName              = "/"+iparam;
FolderProjPath          = BasicPath + FolderName;

nflag0                  = str("%.5d"%cs.Nx)+ "_kyIndex_" + str( "%.5d"%jparam ) + ".dat";

fname0                  = "OccupationCBMath_ncNx_" + nflag0;
fname1                  = "CoherenceMathRealPiNx_" + nflag0;
fname2                  = "CoherenceMathImagPiNx_" + nflag0;

path0                   = FolderProjPath + fname0;
path1                   = FolderProjPath + fname1;
path2                   = FolderProjPath + fname2;


#####################################################
## Creating figure dir.
FigureDir         = '/Figure';
try:
    os.mkdir(BasicPath+FigureDir)#0777);
except OSError as exc:
    if (exc.errno != errno.EEXIST):
        raise exc;
        pass


#####################################################
## Data @################
#x           = []
#f           = open ( path3 , 'r' )
#for line in f:
#    x.append(line);
#f.close()




#####################################################
### LOADING DATA #################
occupationCB_nc     = np.loadtxt( path0 );
coheren_pi          = np.loadtxt( path1 ) + cs.I*np.loadtxt( path2 );



######################################################
## Occupation nc ##
occupationCB_nc = np.transpose( occupationCB_nc );
coheren_pi      = np.transpose( coheren_pi );

print( "\nshape of occup        = ", occupationCB_nc.shape, ";  at ky-index = ", jparam)
print( "\nshape of coherence    = ", coheren_nc.shape, ";  max(cohe-pi)  =  ", coheren_nc.max(), "\n")



#####################################################
## Want to plot occupation nc(kx,ky*,t), need then, kx, ky and t
ky              = -cs.ykmax + cs.dky*jparam;
kx              =  cs.grid( cs.dkx , cs.Nx  );
tme             =  cs.grid( cs.dt , cs.Ntime );



#####################################################
###

for ktime in range(cs.Ntime):        
    
    if (ktime%skiper==0):
        oname = "/nc_it" + str("%.6d"%ktime) + "_iky" +str("%.6d"%jparam) + ".dat";
        qpath = BasicPath + oname;
        
        np.savetxt(
                   qpath
                   ,occupationCB_nc[:,ktime]
                   ,'%.16f'
                   );



#plt.show();


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
