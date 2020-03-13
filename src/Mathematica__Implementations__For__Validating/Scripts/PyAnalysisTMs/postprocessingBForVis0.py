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
#ky              =   cs.grid( cs.dky,  cs.Ny   );
ky = np.zeros(cs.Ny);
for i in range(cs.Ny):
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


t   = laser[:,1]/cs.T0;
set = range(kparam,nparam)

for n in set:
    jparam              = cs.skiper*n
    fname0              = "nc_it"+str("%.6d"%jparam)+".dat"
    path0               = FolderProjPath + fname0;


    #####################################################
    ### LOADING DATA #################
    occupationCB_nc     = np.reshape( np.loadtxt( path0 ), (cs.Ny, cs.Nx) );
    occupationCB_nc     = abs(occupationCB_nc) + 1e-14 ;


    ######################################################
    ## Occupation nc ##
    #occupationCB_nc = np.transpose(occupationCB_nc)
    zmax0   = occupationCB_nc.max()
    zmin0   = occupationCB_nc.min()
    print( "\nMin(nc) = ",zmin0,"Max(nc) = ", zmax0, ";  at time-index = ", jparam)

    zmin0   = 1e-8  #1e-3#1e-14;#-0.055;#
    
    
    #####################################################
    ### Visualization function
    scale   = "log"#"linear"#
    zmax    = 0.442340048601#0.45
    zmin    = zmin0
    fparams = [-1.91, 1.91, -1.1, 1.1, zmin, zmax];
    oparams = [ 0.15,0.1
               ,"kx ","ky"
               , "time = "+str("%.3f"%(jparam*cs.dt/cs.T0))+ " o.c."
               ,"Yes"]

    plt     = cs.occupation_vis(kx,ky,occupationCB_nc,fparams,oparams,scale);



    plt.axes([.152,.1,.25,.2])
    plt.plot( t, laser[:,2], 'r', lw=.35, alpha=0.7 );
    plt.plot( t, laser[:,2], '--r', lw=.15 );

    Ntaux       = int(jparam)+1;
    taux        = np.zeros(Ntaux);
    efaux       = np.zeros(Ntaux);
    y2          = np.zeros(Ntaux);
    
#for i in range(0,Ntaux):
    aix         = np.arange( Ntaux );
    taux        = laser[aix,1]/cs.T0;
    efaux       = laser[aix,2];

    #plt.plot(taux, efaux, 'ro', lw=2 );
    #plt.fill(taux, efaux, 'b', taux, y2, 'b', alpha=0.3)
    plt.fill_between(taux, efaux,y2=0, color=  'g', alpha=0.5)

    #plt.xlabel( 'time', fontsize=8 );
    #plt.tick_params(labelsize=8);
    plt.xlim( .9*t.min(),.90*t.max() );
    plt.ylim( -1.075* laser[:,2].max()
             ,+1.075* laser[:,2].max() );

    plt.xticks([]);
    plt.yticks([]);


    figname = "/"+scale+"Evolution_Occupation_" + str("%.5d"%n) + ".png";
    fileNamePicture     = BasicPath  + FigureDir + figname;


    plt.savefig(fileNamePicture, dpi = 300);
    #shutil.copyfile( fileNamePicture, BasicPath  + figname );




#plt.show();


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
