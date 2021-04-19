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

##################################
### INPUT PARAMETS TO THE SCRIPT
set_of_params   = sys.argv;


iparam          = set_of_params[1]       # give a directory name connected to the path
jparam          = int(set_of_params[2])  # ky-index to be read and analized


################################
##  Basic current path
BasicPath       = os.getcwd();  ##Current directory path




######################################################
################################
#Need to load files
# OccupationCBMath_ncNx_00250_kyIndex_
# CoherenceMathRealPiNx_00250_kyIndex_
FolderName              = "/"+iparam;
FolderProjPath          = BasicPath + FolderName;

nflag0                  = str("%.5d"%cs.Nx)+ "_kyIndex_" + str( "%.5d"%jparam ) + ".bin";
nflag2                  = str("%.5d"%cs.Nx)+ "_kyIndex_" + str( "%.5d"%jparam ) + ".dat"
# file flag for ky-index parameter and *.dat data type
# file flag for ky-index parameter and binary data type

fname0                  = "OccupationCBMath_ncNx_" + nflag0;  #File name for occupation
fname1                  = "CoherenceMathRealPiNx_" + nflag0;
fname2                  = "CoherenceMathImagPiNx_" + nflag2;
fname3                  = "outlaserdata.dat";

path0                   = FolderProjPath + fname0;
path1                   = FolderProjPath + fname1;
path2                   = FolderProjPath + fname2;
path3                   = FolderProjPath + fname3;

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
### LOADING DATA #################
data0            = np.fromfile( path0, dtype = np.dtype("f8") ); #loading occupation


rdata1           = np.fromfile( path1, dtype = np.dtype("f8") )
idata1           = np.loadtxt(  path2  ); #loading imaginary part of coherence
data2            = np.loadtxt(  path3  ); #loading electric field

#idata1           = np.fromfile( path2, dtype = np.dtype("f8") );

######################################################
# Reshaping the occupation and coherence data, data0, rdata1, idata1
occupationCB_nc     = np.reshape( data0,   (cs.Ntime, cs.Nx) );
coheren_pi          = np.reshape( rdata1 , (cs.Ntime, cs.Nx) ) + cs.I*idata1; #+ cs.I*idata1


#Passing the laser field values of data2 to variable efield
efield              = np.zeros( (cs.Ntime,2) );
efield[:,0]         = data2[:,2];




######################################################
## Re-structurating data ##
occupationCB_nc     = np.transpose( occupationCB_nc );  #Transpose of the occupation matrix
coheren_pi          = np.transpose( coheren_pi );       #Transpose of the coherence matrix




############################################################################
#  The current coherence_pi and occupationCB_nc array are 2D, defined by:
#  coherence_pi[ xindex, timeindex]
#  xindex     -> runs on x-direction for kx
#  timeindex  -> runs on time-steps
# The same data structure for occupationCB_nc
######################################


print( "\nshape of occup        = ", occupationCB_nc.shape, ";  at ky-index = ", jparam)
print( "\nshape of coherence    = ", coheren_pi.shape, ";  max(cohe-pi)  =  ", coheren_pi.max(), "\n")



#####################################################
## Calculating momentum axis along x-direction
kx              =  cs.grid( cs.dkx , cs.Nx      );

## Calculating the momentum point ky for the
## parameter index; jparam
ky              = -cs.ykmax + cs.dky*jparam;



######################################################
## Creating time axis
tme             =  cs.grid( cs.dt  , cs.Ntime   );
adipole0        =  np.zeros( (cs.Nx,2) ) + cs.I*np.zeros( (cs.Nx,2) );
agroup_vel_vb   =  np.zeros( (cs.Nx,2) );
agroup_vel_cb   =  np.zeros( (cs.Nx,2) );
aberry_curva_cz  =  np.zeros( cs.Nx );



for i in range(cs.Nx):
    adipole0[i,:]       = cs.xydipoleCV2(   kx[i], ky);
    agroup_vel_cb[i,:]  = cs.xygroup_velC2( kx[i], ky );
    agroup_vel_vb[i,:]  = cs.xygroup_velV2( kx[i], ky );
    aberry_curva_cz[i]  = cs.zBerryCurvaC2( kx[i], ky  );


## Symmetry condition on the dipole
#if ky <= 0:
#    adipole0*=-1.;


######################################################
####
interDip        = np.zeros( (cs.Ntime,2) );
intraJ          = np.zeros( (cs.Ntime,2) );

print ("\n\nCalculation of expectation values along kx as function of time, and at given fixed ky")
print("jparam       = ",jparam,";     ky    =   ", ky)

for ktime in range( cs.Ntime ):
    
    #Calculating inter band contribution
    xtemp1          = -2.*np.real(np.dot( np.conj(adipole0[:,0])
                            ,coheren_pi[:,ktime]   ))*cs.dkx;
    ytemp1          = -2.*np.real(np.dot( np.conj(adipole0[:,1])
                            ,coheren_pi[:,ktime]   ))*cs.dkx;
    
    #Calculating intra band contribution
    xtemp2          = -(  np.dot( agroup_vel_cb[:,0], occupationCB_nc[:,ktime]      )*cs.dkx
                        + np.dot( agroup_vel_vb[:,0], (1.-occupationCB_nc[:,ktime]) )*cs.dkx);
    
    ytemp2          = -(np.dot( agroup_vel_cb[:,1],   occupationCB_nc[:,ktime]      )
                        + np.dot( agroup_vel_vb[:,1], (1. - occupationCB_nc[:,ktime]) ) )*cs.dkx;
    
    #Calculating anomalous velocity
    yanomalousc     =  efield[ktime,0]*aberry_curva_cz
    yanomalousv     = -efield[ktime,0]*aberry_curva_cz
    
    
    ytemp3          = +(np.dot( yanomalousc,occupationCB_nc[:,ktime])*cs.dkx
                       + np.dot(yanomalousv,(1.-occupationCB_nc[:,ktime]) )*cs.dkx)
    

    interDip[ktime,0]   = xtemp1;
    interDip[ktime,1]   = ytemp1;
    intraJ[ktime,0]     = xtemp2;
    intraJ[ktime,1]     = ytemp2 + ytemp3;


oname0      = "/interband_dipole"     + "_iky" +str("%.6d"%jparam) + ".dat";
oname1      = "/intraband_current"    + "_iky" +str("%.6d"%jparam) + ".dat";

qpath0      = BasicPath + oname0;
qpath1      = BasicPath + oname1;

np.savetxt( qpath0, interDip, '%.16f' );
np.savetxt( qpath1, intraJ  , '%.16f' );
#
#
##plt.show();
#
#
print( "\n\n\n/*****************************/" )
print( "Python Program Has finished\n"          )
print( "Hasta la vista BB\n*************\n")
