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


fname0                  = "OccupationCBMath_ncNx_"+str("%.5d"%cs.Nx)+"_kyIndex_" + str( "%.5d"%jparam ) + '.bin';#'.dat';
path0                   = FolderProjPath + fname0;




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
#occupationCB_nc     = np.loadtxt( path0 );
data                = np.fromfile( path0, dtype = np.dtype("f8") );
occupationCB_nc     = np.reshape( data, (cs.Ntime,cs.Nx) );


######################################################
## Occupation nc ##
occupationCB_nc = np.transpose(occupationCB_nc)


print( "\nshape of occup = ", occupationCB_nc.shape, ";  at ky-index = ", jparam, "\n")



#####################################################
## Want to plot occupation nc(kx,ky*,t), need then, kx, ky and t
ky              = -cs.ykmax + cs.dky*jparam;
kx              =  cs.grid( cs.dkx , cs.Nx  );
tme             =  cs.grid( cs.dt , cs.Ntime );




#jparam          = jparam+1
#k=0

#####################################################
## Sypers
zmin0 = occupationCB_nc.min()
zmax0 = occupationCB_nc.max()


print ("\nSkyper    =   ",   cs.skiper ,";     Ntime'  =   ", cs.skiper*cs.Nsnap,";      Real Ntime  = ",cs.Ntime )
print ("ky-index    =  ",jparam,";           max(n_c)    =   ",zmax0,";            Min of nc = ",zmin0 )




#####################################################
## OutPut for Occupation dynamics
for ktime in range(cs.Ntime):
    if (ktime%cs.skiper==0):
        oname = "/nc_it" + str("%.6d"%ktime) + "_iky" +str("%.6d"%jparam) + ".dat";
        qpath = BasicPath + oname;
        
        np.savetxt(
                   qpath
                   ,occupationCB_nc[:,ktime]
                   ,'%.16f'
                   );




#####################################################
### Visualization function
zmax    =   +0.14;
zmin    =   zmin0;

scale       =   "linear";
axesparam   =   [-8,    8,   -1.91,  1.91,   zmin,  zmax];
oparams     =   [0.1,   0.1,    "time ",    "kx",   "ky = " +str("%.3f"%ky), "Yes"];



#####################################################
## Occupations ....
#plt     = cs.occupation_vis( tme/cs.T0
#                            ,kx
#                            ,occupationCB_nc
#                            ,axesparam
#                            ,oparams
#                            ,scale
#                            );
#
#
#
#figname = "/" + scale + "Occupation_Density_ky_" + str("%.5d"%jparam) + ".png";
#fileNamePicture     = BasicPath  + FigureDir + figname;
#
#
#plt.savefig(fileNamePicture, dpi = 300);
#shutil.copyfile( fileNamePicture, BasicPath  + figname );
#plt.show();


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
