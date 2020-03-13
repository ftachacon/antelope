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
#Nsnap   = 160; #snaptshot per opt. cycle
#snaper  = 16;#32

#####################################################
## Want to plot occupation nc(kx,ky*,t), need then, kx, ky and t
kx              =   cs.grid( cs.dkx , cs.Nx   );
ky              =   np.zeros(cs.Ny+1);
for i in range(cs.Ny+1):
    ky[i]   = -cs.dky*(cs.Ny-1.)/2. + cs.dky*i;

tme             =   cs.grid( cs.dt , cs.Ntime );


################################
### INPUT PARAMETS
set_of_params   = sys.argv;
iparam          = set_of_params[1]

kparam          = int(set_of_params[2]);
nparam          = int(set_of_params[3]);


################################
##  Basic current path
BasicPath       = os.getcwd();



################################
#Need to load files
# OccupationCBMath_ncNx_00250_kyIndex_
# CoherenceMathRealPiNx_00250_kyIndex_
FolderName              = "/"+iparam;
FolderProjPath          = BasicPath + FolderName;



flname              = "outlaserdata.dat";
pathl               = FolderProjPath + flname;
laser               = np.loadtxt( pathl );



#####################################################
## Creating figure dir.
FigureDir               = '/Figure';
try:
    os.mkdir(BasicPath+FigureDir)#0777);
except OSError as exc:
    if (exc.errno != errno.EEXIST):
        raise exc;
        pass

inter_dip_ky = np.zeros( (cs.Ntime-1,2) )
intra_jc_ky  = np.zeros( (cs.Ntime-1,2) )
set          = range(kparam,nparam)

for n in set:
    
    jparam      = n;
    
    fname0      = "/interband_dipole"     + "_iky" +str("%.6d"%jparam) + ".dat";
    fname1      = "/intraband_current"    + "_iky" +str("%.6d"%jparam) + ".dat";
    path0       = FolderProjPath + fname0;
    path1       = FolderProjPath + fname1;
    
    
    #####################################################
    ### LOADING DATA #################
    inter_dip =  np.loadtxt( path0 );
    intra_cur =  np.loadtxt( path1 );
    
    
    #np.reshape( np.loadtxt( path0 ), (cs.Ny+1, cs.Nx) );
    inter_dip_ky[:,0]+= inter_dip[:,0]*cs.dky;
    inter_dip_ky[:,1]+= inter_dip[:,1]*cs.dky;

    intra_jc_ky[:,0]+=  intra_cur[:,0]*cs.dky;
    intra_jc_ky[:,1]+=  intra_cur[:,1]*cs.dky;

#plt.savefig(fileNamePicture, dpi = 300);
#shutil.copyfile( fileNamePicture, BasicPath  + figname );


oname0      = "/interband_dipole_full_evol.dat";
oname1      = "/intraband_current_full_evol.dat";

qpath0      = BasicPath + oname0;
qpath1      = BasicPath + oname1;

np.savetxt( qpath0, inter_dip_ky, '%.16f' );
np.savetxt( qpath1, intra_jc_ky, '%.16f' );




width   = 8
hight   = width/1.62
ax01     = 0.1
ax02     = 0.05
fig     = plt.figure(figsize=(width,hight));
ax1     = fig.add_axes([ax01, ax02, 0.85, 0.85])

#plt.plot(laser[0:-2,1],laser[0:-2,2])
t   = laser[0:-1,1]/cs.T0
plt.plot(t,laser[0:-1,2],"r",label='$E_{L}$')
plt.plot(t,inter_dip_ky[:,0],"b",label='$J_{x,er}$')
plt.plot(t,inter_dip_ky[:,1],"g",label='$J_{y,er}$')

plt.legend();


plt.title("inter-band, x and y compontent")
print("\nlasershape",laser.shape, '; time shape = ', t.shape )



fig     = plt.figure(figsize=(width,hight));
ax1     = fig.add_axes([ax01, ax02, 0.85, 0.85])
plt.plot( t,laser[0:-1,2],"r",label='$E_{L}$')
plt.plot( t,intra_jc_ky[:,0],"b",label='$J_{x,ra}$')
plt.plot( t,intra_jc_ky[:,1],"g",label='$J_{y,ra}$')

plt.legend();


plt.title("intra-band, x and y compontent")

plt.show();


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
