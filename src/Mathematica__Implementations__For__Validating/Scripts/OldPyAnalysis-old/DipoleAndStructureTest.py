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



for n in range(3):
    print ("\n\na", n, " = (", cs.vec_a[0][n],",  ", cs.vec_a[1][n]," )")
    print ("\nb", n, " = (", cs.vec_b[0][n],",  ", cs.vec_b[1][n]," )")




################################
### INPUT PARAMETS
set_of_params   = sys.argv;
iparam          = set_of_params[1]



################################
##  Basic current path
BasicPath       = os.getcwd();





################################
#Need to load files
# OccupationCBMath_ncNx_00250_kyIndex_
# CoherenceMathRealPiNx_00250_kyIndex_
FolderName              = "/"+iparam;
FolderProjPath          = BasicPath + FolderName;

#fname0                  = "OccupationCBMath_ncNx_00250_kyIndex_" + str( "%.5d"%jparam ) + '.dat';
#path0                   = FolderProjPath + fname0;


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



#####################################################
## Want to plot occupation nc(kx,ky*,t), need then, kx, ky and t
kx              =  cs.grid( cs.dkx , cs.Nx  );
ky              =  cs.grid( cs.dky , cs.Ny  );#-cs.ykmax + cs.dky*jparam;

kxtest = -1.3
kytest = 0.

chern_no0       = cs.Chern_Number( kx, ky );
print( "\nCNo = ",chern_no0, "\n\n" );

i = 0; j = 0;



ec0     = cs.Ec2( kxtest, kytest );
ev0     = cs.Ev2( kxtest, kytest );

d0      = cs.xydipoleCV2( kxtest, kytest );
bcon0   = cs.xBerryConnectionC2( kxtest, kytest);
bcon0   = cs.yBerryConnectionC2( kxtest, kytest );
vg_v0   = cs.xygroup_velV2( kxtest, kytest );
vg_c0   = cs.xygroup_velC2( kxtest, kytest );
bcurva  = cs.zBerryCurvaC2( kxtest, kytest  );


print ("\n\nEnergy Dispersion, Ec = ", ec0, ";  at kx =  ", kxtest,";  at ky =  ", kytest )
print ("Energy Dispersion, Ev = ", ev0, ";  at kx =  ", kxtest,";  at ky =  ", kytest, "\n" )
print ("\nDipoleCV = ", d0, ";  at kx =  ", kxtest, ";  at ky =  ", kytest, "\n" )
print ("\nxBConnec = ", bcon0, ";  at kx =  ", kxtest,";  at ky =  ", kytest )
print ("\nyBConnec = ", bcon0, ";  at kx =  ", kxtest,";  at ky =  ", kytest )
print ("\nvg_v = ", vg_v0, ";  at kx =  ", kxtest,";  at ky =  ", kytest, "\n" )
print ("\nvg_c = ", vg_c0, ";  at kx =  ", kxtest,";  at ky =  ", kytest, "\n" )


Bcurva = np.zeros( (cs.Ny,cs.Nx) )
for j in range(cs.Ny):
    for i in range(cs.Nx):
        temp=cs.xydipoleCV2(kx[i],ky[j])
        Bcurva[j,i] = np.real(temp[1]) #zBerryCurvaC2( kx[i],ky[j] )



#####################################################
### Visualization function
#          xmin     xmax    ymin   ymax    zmin      zmax
fparams = [-1.91,    1.91,   -1.1,  1.1,   Bcurva.min(),  Bcurva.max()];
oparams = [0.15,0.1,"kx","ky", "y-dir dipole re "]


plt     = cs.occupation_vis(kx*cs.T0,ky,Bcurva,fparams,oparams);


#figname = "/BCurva_" + ".png";
figname = "/yDipoleReC_" + ".png";
fileNamePicture     = BasicPath  + FigureDir + figname;


plt.savefig(fileNamePicture, dpi = 300);
#shutil.copyfile( fileNamePicture, BasicPath  + figname );



plt.show();


print( "\n\n\n/*****************************/")
print( "Python Program Has finished\n")
print( "Hasta la vista BB\n*************\n")
