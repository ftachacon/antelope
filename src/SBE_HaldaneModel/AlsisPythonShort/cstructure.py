#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:42:05 2019

@author: achacon
"""


#//  solidstructure.h
#//  Created by Alexis Agustín  Chacón Salazar on 3/19/19.
#// Structural information defined as follow:

#/***********************************************
 
# (1) Energy dispersions,
# (2) Dipole matrix elements,
# (3) Berry connection and
# (4) Berry curvature,
# (5) Group velocities,
# (6) Anomalous velocities
# (7) Band gap
# (8) Chern No.
 
#**********************************/

#ifndef solidstructure_h
#define solidstructure_h

#include <stdlib.h>
#include <string>
#include "mkl.h"
#include "constant.h"
#include "momaxis.h"
#using namespace std;

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from matplotlib.colors import LogNorm


import numpy as np
import os 
import sys 
import shutil 
import cmath

from scipy.fftpack import fft, ifft
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.interpolate import interp2d

pi      = np.pi;
I       = 1j;


###############################
Nvects      = 3;
Ndim        = 2;
Nsnap       = 160*2; #snaptshot per opt. cycle


###############################
#Lattice constant
a0          = 1./0.529;
t1          = 0.075;
t2          = t1/3.;
M0          = 2.54*t2;
phi0        = 1.16;
eps_phi0    = 0.5e-5;



###############################
## MOMENTUM AXIS
Nx      = 400;      #801; #401; #301;#401;#301;#101;
Ny      = 100;      #100; #101; #;

xkmax   = 1.918999727359802;
ykmax   = 1.107935009166000;

dkx     = 2.*xkmax/( Nx - 1. );
dky     = 2.*ykmax/( Ny - 1. );

################################
## TIME STEPS
dt      = 2.;#1.20; # 1.61; #1.;
#1.1219973762#1.401970981235;
#1.869995627136782;#1.4024967203525864
#2.804993440705173;#1.401970981235;#     #1.401970981234623283439

Ntime   = 7180; #11968#10771#8616#10053#10053#11489#7180#8975; #14361#10771#15387;
#10258;#15387#14361;#12799#10240#11520; #9600;
#10240#5120#4480#3840#5120##2560;#640#3200;      #3585#5121;#       #2560;             

skiper  = int(Ntime/Nsnap);

w0      = 0.014;
T0      = 2.*pi/w0;
E0      = 0.006;                #0.007

#############################
#bases_vectors()
# First vector base: position atom
vec_a       = np.zeros( (2,Nvects) )

vec_a[0,0]  = 0.;
vec_a[1,0]  = a0;

vec_a[0,1]  = -np.sqrt(3.)/2.*a0;
vec_a[1,1]  = -1./2.*a0;

vec_a[0,2]  =  np.sqrt(3.)/2.*a0;
vec_a[1,2]  = -1./2.*a0;

#b-vectors of the hexagonal lattice
vec_b       = np.zeros( (2,Nvects) )
vec_b[0,0]  = np.sqrt( 3. )*a0;
vec_b[1,0]  = 0.;

vec_b[0,1]  = -np.sqrt(3.)/2.*a0;
vec_b[1,1]  = +3./2.*a0;

vec_b[0,2]  = -np.sqrt(3.)/2.*a0;
vec_b[1,2]  = -3./2.*a0;


#//K and K' points
K           = np.zeros(2);
K[0]        = +4.*pi/np.sqrt(3.)/3./a0 ;
K[1]        = 0. ;

Kprime      = np.zeros(2)
Kprime[0]   = -4.*pi/np.sqrt(3.)/3./a0 ;
Kprime[1]   = 0. ;



def grid(dk,N):
    k = np.zeros(N);
    kmax = (N-1)*dk/2.;
    
    for n in range(0,N):
        k[n] = -kmax + n*dk;
    
    return k

#######################################################
############### Visualization tool ################
def occupation_vis(time,kx,occupationCB_nc,fparams,oparams,scale):

    xmin    = fparams[0]
    xmax    = fparams[1]
    
    ymin    = fparams[2]
    ymax    = fparams[3]
    
    ax1 = oparams[0]
    ax2 = oparams[1]
    
    xlab = oparams[2]
    ylab = oparams[3]
    
    atitle = oparams[4]
    
    
    X,Y     = np.meshgrid(time,kx);
    
    
    #f       = interp2d( X, Y, occupationCB_nc, kind='cubic')
    #xnew    = np.arange(kx[0], kx[-1], 1*Nx)
    #ynew    = np.arange(ky[0], ky[-1], 2*Ny)
    #data1   = f(xnew,ynew)
    #Xn, Yn  = np.meshgrid(xnew, ynew)
    
    #print("Xnshape,", Xn.shape)
    
    nmin    = fparams[4]#occupationCB_nc.min(); #np.log10( occupationCB_nc.min())
    nmax    = fparams[5]#occupationCB_nc.max(); #np.log10( occupationCB_nc.max() );


    levels = MaxNLocator(nbins=18).tick_values( nmin, nmax );


    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.

    colorMap0='bwr'#'coolwarm'#'gist_rainbow'#'RdBu'
    #'seismic'#'gray_r'#'RdGy'#'PRGn'#'bwr'#'gnuplot'#'gnuplot2'#'CMRmap_r'#'coolwarm'#'CMRmap_r'
    cmap = plt.get_cmap(colorMap0)


    norm = BoundaryNorm(levels
                    ,ncolors=cmap.N
                    ,clip=True
                    );


    #####################################################
    ## Visualization
    fontz  = 30;

    width = 12
    hight = width/1.62

    fig = plt.figure(figsize=(width,hight));
    ax1 = fig.add_axes([ax1, ax2, 0.85, 0.85])


#cf = plt.contourf( T
#                  ,Q
#                  ,occupationCB_nc
#                  ,levels=levels
#                  ,cmap=cmap );
#
#pcolormesh
    if scale=="linear":
        cf=plt.pcolormesh( X
                      ,Y
                      ,occupationCB_nc
                      ,vmin=nmin
                      ,vmax=nmax #,levels=levels
                      ,cmap=cmap );
        
                      
    if scale=="log":
            cf=plt.pcolormesh( X
                                ,Y
                                ,occupationCB_nc
                                ,norm = LogNorm( vmin=nmin
                                ,vmax=nmax )
                                ,cmap=cmap );
 
 
    if (oparams[5]=="Yes"):
       cb = plt.colorbar(cf);#mappable
    #else:
    #cb = plt.colorbar(cf);#mappable
    #cb.ax.set_ticks([]);
    #cb.ax.tick_params([]);
    #cb.ax.tick_params(labelsize=20);
    #ax1.set_title('contourf with levels')
    #plt.title('contourf with levels',fontsize=20)

    plt.title(atitle,fontsize=21 )

    plt.xlabel(xlab, fontsize=21);
    plt.ylabel(ylab, fontsize=21);

    plt.tick_params(labelsize=21);

#   plt.xticks(np.arange(-2,2,0.5))
#   plt.yticks(np.arange(-1.,1.5,0.5))

    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])

    return plt

    
###########################
# Dipole matrix element
# fx -> pi_y(kx)
def inter_dipole(kx,ky,fx):
    integ = np.zeros(2)+I*np.zeros(2);
    dipoleVec
    for n in range(Nx):
        dip_cv  =  xydipoleCV2( kx[n], ky )
        integ[0]+= fx[n]*dip_cv[0]
        integ[1]+= fx[n]*dip_cv[1]
    integ[0]*= -dkx;
    integ[1]*= -dkx;
    return integ


#####################################
#####################################
#######  group velocity average #####
def jvg(kx,ky,nc):
# jvg[0] = nc*vg_c[0] + nv*vg_v[0]
#
    Nx  = len(kx);
    dkx = kx[1]-kx[0];
    jv  = np.zeros(2);
    
    for n in range(Nx):
        gvelv = xygroup_velV2( kx[n], ky )
        gvelc = xygroup_velC2( kx[n], ky )
        
        jv[0]+= nc[n]*gvelc[0] + (1.-nc[n])*gvelv[0]
        jv[1]+= nc[n]*gvelc[1] + (1.-nc[n])*gvelv[1]

    jv*=    dkx;
    return jv;





#####################################s
#####################################
#### anomalous velocity average  ####
def jva(kx, ky, nc, efield ):
    # jvg[0] = nc*va_c[0] + nv*va_v[0]
    #
    Nx  = len(kx);
    dkx = kx[1]-kx[0];
    ja  = np.zeros(2);
    
    for n in range(Nx):
        
        zcurvac0    =  zBerryCurvaC2( kx[n], ky  );
        zcurvav0    = -zcurvac0;
        
        janalc0     = xyanomalous_velC( efield, zcurvac0 );
        janalv0     = xyanomalous_velV( efield, zcurvav0 );
        
        ja[0]+= (nc[n]*janalc0[0] + (1.-nc[n])*janalv0[0])
        ja[1]+= (nc[n]*janalc0[1] + (1.-nc[n])*janalv0[1])
    
    ja*= -dkx;
    return ja;




#####################################
#####################################
## Intra band current calculation  ##
def j_ra(kx,ky,nc,efield):

    jgroupv0 = jvg(kx,ky,nc);
    janomal0 = jva(kx,ky,nc,efield);

    return -(jgroupv0+janomal0);


class Ham:
    
    def __init__(self):
        
        self.b0res0         = 0.;
        self.b1res0         = 0.;
        self.b2res0         = 0.;
        self.b3res0         = 0.;

    
        self.xb0der0        = 0.;
        self.xb1der0        = 0.;
        self.xb2der0        = 0.;
        self.xb3der0        = 0.;
    
        self.yb0der0        = 0.;
        self.yb1der0        = 0.;
        self.yb2der0        = 0.;
        self.yb3der0        = 0.;
    
        self.bnorm0             = 0.;
        self.xbnorm0            = 0.;
        self.ybnorm0            = 0.;
    
        self.phibloch0          = 0.;
        self.thetabloch0        = 0.;
    
        self.xgradbnorm0        = 0.;
        self.ygradbnorm0        = 0.;

        self.xgrad_phibloch0    = 0.;
        self.ygrad_phibloch0    = 0.;

        self.xgrad_thetabloch0  = 0.;
        self.ygrad_thetabloch0  = 0.;
    
        self.xgrad_phibloch02   = 0.;
        self.ygrad_phibloch02   = 0.;
    
        self.xgrad_thetabloch02 = 0.;
        self.ygrad_thetabloch02 = 0.;
    
        self.ec0                = 0.;
        self.ev0                = 0.;
        self.eg0                = 0.;
    
        self.xberryconnCB0      = 0.;
        self.xberryconnVB0      = 0.;

        self.yberryconnCB0      = 0.;
        self.yberryconnVB0      = 0.;

        self.xchig0             = 0.;
        self.ychig0             = 0.;
    
        self.xdip_cv0           = 0.*I;
        self.ydip_cv0           = 0.*I;
    
        self.xgroup_c0          = 0.;
        self.ygroup_c0          = 0.;

        self.xgroup_v0          = 0.;
        self.ygroup_v0          = 0.;

        self.zberry_curvaCB0    = 0.;
        self.zberry_curvaVB0    = 0.;


    ############################################
    ############################################
    # Components of the Haldane model Hamiltonian
    def B0( self, kx, ky ):
        self.b0res0 = 0.;
        for n in range(Nvects):
            self.b0res0+= np.cos( kx*vec_b[0,n] + ky*vec_b[1,n] );
        self.b0res0*= 2. * t2 * np.cos( phi0 );
        self.b0res0;

    ######################################
    def B1( self, kx, ky ):
        self.b1res0 = 0.;
        for n in range( Nvects ):
            self.b1res0+= np.cos( kx*vec_a[0,n] + ky*vec_a[1,n] );
        self.b1res0*= t1;
        self.b1res0;

    def B2( self, kx, ky ):
        self.b2res0 = 0.;
        for n in range( Nvects ):
            self.b2res0+= np.sin( kx*vec_a[0,n] + ky*vec_a[1,n] );
        self.b2res0*= t1;
        self.b2res0;

    def B3( self, kx, ky ):
        self.b3res0 = 0.;
        for n in range( Nvects ):
            self.b3res0+= np.sin( kx*vec_b[0,n] + ky*vec_b[1,n] );
        
        self.b3res0*= -2.0 * t2 * np.sin( phi0 );
        self.b3res0+= M0;
        self.b3res0;

    def BNorm( self ):
        self.bnorm0 = np.sqrt(self.b1res0 * self.b1res0+self.b2res0 * self.b2res0+self.b3res0 * self.b3res0);
    
    def phiBloch(  self ):
        self.phibloch0 = np.arctan2( self.b2res0, self.b1res0 );
    
    def thetaBloch( self ):
        self.thetabloch0= np.arccos( self.b3res0 / self.bnorm0 );
    
    def xgradB0( self, kx, ky ):
        self.xb0der0 = 0.;
        for n in range( Nvects ):
            self.xb0der0+=  vec_b[0,n]*np.sin( vec_b[0,n]*kx + vec_b[1,n]*ky );
        self.xb0der0*= -2.*t2*np.cos( phi0 );

    def xgradB1( self, kx, ky ):
        self.xb1der0 = 0.;
        for n in range(0, Nvects ):
            self.xb1der0+=  vec_a[0,n]*np.sin( vec_a[0,n]*kx + vec_a[1,n]*ky );
        self.xb1der0*= -t1;

    def xgradB2( self, kx, ky ):
        self.xb2der0 = 0.;
        for n in range( Nvects ):
            self.xb2der0+=  vec_a[0,n]*np.cos( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );
        self.xb2der0*= t1;

    def xgradB3( self, kx, ky ):
        self.xb3der0 = 0.;
        for n in range( Nvects ):
            self.xb3der0+=  vec_b[0,n]*np.cos( vec_b[0,n]*(kx) + vec_b[1,n]*(ky) );
        self.xb3der0*= -2.*t2*np.sin( phi0 );

    def ygradB0( self, kx, ky ):
        self.yb0der0 = 0.;
        for n in range( Nvects ):
            self.yb0der0+=  vec_b[1,n]*np.sin( vec_b[0,n]*kx + vec_b[1,n]*ky );
        self.yb0der0*= -2.*t2*np.cos( phi0 );

    def ygradB1(self, kx, ky ):
        self.yb1der0 = 0.;
        for n in range( Nvects ):
            self.yb1der0+=  vec_a[1,n]*np.sin( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );
        self.yb1der0*= -t1;

    def ygradB2( self, kx, ky ):
        self.yb2der0 = 0;
        for n in range( Nvects ):
            self.yb2der0+=  vec_a[1,n]*np.cos( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );
        self.yb2der0*= t1;

    def ygradB3( self, kx, ky ):
        self.yb3der0 = 0.;
        for n in range( Nvects ):
            self.yb3der0+=  vec_b[1,n]*np.cos( vec_b[0,n]*(kx) + vec_b[1,n]*(ky) );
        self.yb3der0*= -2.*t2*np.sin( phi0 );
    
    def xgradBNorm( self ):
        self.xbnorm0 = ( self.xb1der0*self.b1res0 + self.xb2der0*self.b2res0 + self.xb3der0*self.b3res0 )/self.bnorm0;


    def ygradBNorm( self  ):
        self.ybnorm0 = ( self.yb1der0*self.b1res0 + self.yb2der0*self.b2res0 + self.yb3der0*self.b3res0 )/self.bnorm0;

    
    def xgrad_phiBloch2( self, eps1 ):
            self.xgrad_phibloch02 = ( self.b1res0 * self.xb2der0 - self.b2res0 * self.xb1der0 )/( self.b1res0 * self.b1res0 + self.b2res0 * self.b2res0 + eps1*eps1 );
            self.xgrad_phibloch0 = ( self.b1res0 * self.xb2der0 - self.b2res0 * self.xb1der0 )/( self.b1res0 * self.b1res0 + self.b2res0 * self.b2res0 + 1.e-40 );


    def xgrad_thetaBloch2( self, eps1 ):
        temp =  self.bnorm0*self.bnorm0 - self.b3res0*self.b3res0;
        
        self.xgrad_thetabloch02 = -( self.bnorm0*self.xb3der0 - self.xbnorm0*self.b3res0 )/(np.sqrt( temp ) + eps1)/(self.bnorm0 + eps1);
        
        
        self.xgrad_thetabloch0 = -( self.bnorm0 * self.xb3der0 -  self.xbnorm0 * self.b3res0 ) / (np.sqrt( self.bnorm0*self.bnorm0 - self.b3res0*self.b3res0 ) + 1.e-20)/ (self.bnorm0 + 1.e-20);
    

    def ygrad_phiBloch2(  self, eps2  ):
        self.ygrad_phibloch02= ( self.b1res0 * self.yb2der0 - self.b2res0 * self.yb1der0 )/( self.b1res0 * self.b1res0 + self.b2res0 * self.b2res0 + eps2  );
        self.ygrad_phibloch0= ( self.b1res0 * self.yb2der0 - self.b2res0 * self.yb1der0 )/( self.b1res0 * self.b1res0 + self.b2res0 * self.b2res0 + 1.e-20  );

    def ygrad_thetaBloch2( self, eps2 ):
        temp = self.bnorm0 * self.bnorm0 - self.b3res0 * self.b3res0
        self.ygrad_thetabloch02 = -( self.bnorm0 * self.yb3der0 -  self.ybnorm0 * self.b3res0 )/( np.sqrt( temp  ) + eps2)/(self.bnorm0  + eps2);
        self.ygrad_thetabloch0 = -( self.bnorm0 * self.yb3der0 -  self.ybnorm0 * self.b3res0 )/( np.sqrt( temp  ) +1.e-20)/(self.bnorm0 +1.e-20);

    def bset(self, kx, ky):
        self.B0(kx,ky);
        self.B1(kx,ky);
        self.B2(kx,ky);
        self.B3(kx,ky);
        self.BNorm();
        self.phiBloch();
        self.thetaBloch();
    
    def bsetgrads( self, kx, ky, xeps1, yeps1, xeps2, yeps2):
        self.xgradB0(kx,ky);
        self.xgradB1(kx,ky);
        self.xgradB2(kx,ky);
        self.xgradB3(kx,ky);
        self.xgradBNorm(  );
        
        self.ygradB0(kx,ky);
        self.ygradB1(kx,ky);
        self.ygradB2(kx,ky);
        self.ygradB3(kx,ky);
        self.ygradBNorm(  );

        self.xgrad_phiBloch2( xeps1 );
        self.ygrad_phiBloch2( yeps1 );

        self.xgrad_thetaBloch2( xeps2 );
        self.ygrad_thetaBloch2( yeps2 );

    ####################################################
    ####################################################
    # ENERGY DISPERSIONS CONDUCTION AND VALENCE BANDS
    def Ec( self ):
        self.ec0 = self.b0res0 + self.bnorm0;
    def Ev( self ):
        self.ev0 = self.b0res0 - self.bnorm0;
    def Eg( self ):
        self.eg0 = 2.*self.bnorm0;

    ####################################################
    #####################################################
    ## Conduction Band Berry Connection
    def xBerryConnectC( self, eps ):
        self.xberryconnCB0 = self.b3res0/2./( self.bnorm0 + eps)*self.xgrad_phibloch02;

    def xBerryConnectV( self, eps ):
        self.xberryconnVB0 = -self.b3res0/2./( self.bnorm0 + eps)*self.xgrad_phibloch02;

    def yBerryConnectC( self, eps ):
        self.yberryconnectCB0 =self.b3res0/2./( self.bnorm0 + eps)*self.ygrad_phibloch02;

    def yBerryConnectV( self, eps ):
        self.yberryconnectVB0 = -self.b3res0/2./( self.bnorm0 + eps)*self.ygrad_phibloch02;


    ####################################################
    #//Berry Connection difference
    def xyChig( self ):
        self.xchig0  = 2.0 * self.xberryconnCB0;
        self.ychig0  = 2.0 * self.yberryconnCB0;


    ####################################################
    #//Dipole Matrix Element from Val. to Cond. Band
    ##
    def xydipoleCV( self, kx, ky):
    
        if (ky<0):
            self.xgrad_phibloch02*=-1.;
            self.ygrad_thetabloch02*=-1.;

        self.xdip_cv0 =  np.sin( self.thetabloch0 ) * self.xgrad_phibloch02/2. + I * self.xgrad_thetabloch02/2.;
        self.ydip_cv0 =  np.sin( self.thetabloch0 ) * self.ygrad_phibloch02/2. + I * self.ygrad_thetabloch02/2.;
        
        #self.xdip_cv0 = 1 + I
        #self.ydip_cv0 = 1 + I


    ####################################################
    ####################################################
    # GROUP VELOCITIES
    # conduction band group vel.
    def xygroup_velC( self ):
        self.xgroup_c0 =  self.xb0der0 + self.xbnorm0 ;
        self.ygroup_c0 =  self.yb0der0 + self.ybnorm0 ;

    def xygroup_velV( self ):
        self.xgroup_v0 =  self.xb0der0 - self.xbnorm0 ;
        self.ygroup_v0 =  self.yb0der0 - self.ybnorm0 ;


####################################################f
##  BERRY CURVATUREs Conduction band
    def zBerryCurvaC( self ):
        self.zberry_curvaCB0 = - 1./2. * np.sin( self.thetabloch0  ) * ( self.xgrad_thetabloch0 * self.ygrad_phibloch0 - self.xgrad_phibloch0 * self.ygrad_thetabloch0 );

    def zBerryCurvaV( self ):
        self.zberry_curvaVB0 = + 1./2. * np.sin( self.thetabloch0  ) * ( self.xgrad_thetabloch0 * self.ygrad_phibloch0 - self.xgrad_phibloch0 * self.ygrad_thetabloch0 );


    ####################################################
    # Chern No. calculation
    def Chern_Number(self,kx,ky):
        Nx  = len(kx);
        Ny  = len(ky);
        dkx = kx[1]-kx[0];
        dky = ky[1]-ky[0];
        xeps1, yeps1, xeps2, yeps2 = 1.e-16,1.e-16,1.e-16,1.e-16;
        chernNo = 0.;
        for j in range(Ny):
            for i in range(Nx):
                self.bset( kx[i],ky[j] );
                self.bsetgrads( kx[i], ky[j], xeps1, yeps1, xeps2, yeps2 );
                self.zBerryCurvaC( );
                zbcuva = self.zberry_curvaCB0;
                chernNo+= zbcuva;

        return chernNo*dkx*dky/4./pi


#
def xgrad_phiBloch( b1res0, b2res0, xb1der0, xb2der0 ):
    return ( b1res0 * xb2der0 - b2res0 * xb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0 );

def ygrad_rphiBloch( b1res0, b2res0, yb1der0, yb2der0 ):
    return ( b1res0 * yb2der0 - b2res0 * yb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0 + eps_phi0 );
#
def ygrad_phiBloch( b1res0,b2res0,yb1der0, yb2der0  ):
    return ( b1res0 * yb2der0 - b2res0 * yb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0  );

#
#
def xgrad_rthetaBloch( b3res0, bnorm0, xb3der0, xbnorm0 ):
    return ( bnorm0 * xb3der0 -  xbnorm0 * b3res0 )/( np.sqrt( bnorm0 * bnorm0 - b3res0 * b3res0  ) + eps_phi0 )/( bnorm0 + eps_phi0 );

#

def ygrad_rthetaBloch( b3res0, bnorm0, yb3der0, ybnorm0 ):
    return -( bnorm0 * yb3der0 -  ybnorm0 * b3res0 )/( np.sqrt( bnorm0 * bnorm0 - b3res0 * b3res0  ) + eps_phi0 )/( bnorm0 + eps_phi0 );



####################################################
####################################################
# ANOMALOUS VELOCITIES
def xyanomalous_velV( efield, zBCurva_v0  ):
    xanomalous_v0   =  - efield[1] * zBCurva_v0;
    yanomalous_v0   =    efield[0] * zBCurva_v0;

    return [xanomalous_v0,yanomalous_v0]




def xyanomalous_velC( efield, zBCurva_c0 ):
    xanomalous_c0   =  - efield[1] * zBCurva_c0;
    yanomalous_c0   =    efield[0] * zBCurva_c0;

    return [xanomalous_c0,yanomalous_c0]


#
#void solidstructure::energy_vc_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
#{
#
#
#    for ( l = 0; l < g->N[2]/skip2; l++ )
#        for ( j = 0; j < g->N[1]/skip1; j++ )
#        {
#            fprintf( output, "\n");
#            for ( i = 0; i < g->N[0]/skip0; i++ )
#            {
#
#                set_of_B( &g->k[0][i*skip0]
#                         ,&g->k[1][j*skip1] );
#
#                Ev();
#                Ec();
#
#                fprintf( output, "%e %e %.16e %.16e %.16e \n"
#                        ,g->k[0][i*skip0]
#                        ,g->k[1][j*skip1]
#                        ,energyv0
#                        ,energyc0
#                        ,energyc0 - energyv0
#                        );
#
#            }
#
#        }
#
#
#}




#
#void solidstructure::dipole_cv_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
#{
#
#    if (flag==0)
#    {
#
#        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
#        exit (EXIT_FAILURE);
#
#    }
#
#    for ( l = 0; l < g->N[2]/skip2; l++ )
#    for ( j = 0; j < g->N[1]/skip1; j++ )
#    {
#
#        fprintf( output, "\n");
#
#        for ( i = 0; i < g->N[0]/skip0; i++ )
#        {
#
#
#            /*set_of_B( &g->k[0][i*skip0]
#             ,&g->k[1][j*skip1]
#             );
#
#             set_of_Bgrad( &g->k[0][i*skip0]
#             ,&g->k[1][j*skip1]
#             );
#
#             dipoleCV();*/
#
#
#            ai = i*skip0;
#            aj = j*skip1;
#            al = l*skip2;
#
#
#            fprintf( output, "%e %e %.16e %.16e %.16e %.16e \n"
#                    ,g->k[0][i*skip0]
#                    ,g->k[1][j*skip1]
#                    ,real( dipole_cv[0][g->index( &ai, &aj, &al ) ] )
#                    ,imag( dipole_cv[0][g->index( &ai, &aj, &al ) ] )
#                    ,real( dipole_cv[1][g->index( &ai, &aj, &al ) ] )
#                    ,imag( dipole_cv[1][g->index( &ai, &aj, &al ) ] )
#                    );
#
#        }
#
#    }
#
#    fflush(output);
#}



#void solidstructure::connection_c_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
#{
#
#
#    if (flag==0)
#    {
#
#        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
#        exit (EXIT_FAILURE);
#
#    }
#
#    for ( l = 0; l < g->N[2]/skip2; l++ )
#    for ( j = 0; j < g->N[1]/skip1; j++ )
#    {
#
#        fprintf( output, "\n");
#
#        for ( i = 0; i < g->N[0]/skip0; i++ )
#        {
#
#
#
#            ai = i*skip0;
#            aj = j*skip1;
#            al = l*skip2;
#
#            fprintf( output, "%e %e %e %e \n"
#                    ,g->k[0][i*skip0]
#                    ,g->k[1][j*skip1]
#                    ,chi_gap[0][  g->index( &ai, &aj, &al ) ]/2.
#                    ,chi_gap[1][  g->index( &ai, &aj, &al ) ]/2.
#                    );
#
#        }
#
#    }
#
#    fflush(output);
#}


#void solidstructure::curvature_c_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
#{
#
#    if (flag==0)
#    {
#
#        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
#        exit (EXIT_FAILURE);
#
#    }
#
#    for ( l = 0; l < g->N[2]/skip2; l++ )
#    for ( j = 0; j < g->N[1]/skip1; j++ )
#    {
#
#        fprintf( output, "\n");
#
#        for ( i = 0; i < g->N[0]/skip0; i++ )
#        {
#
#
#            ai = i*skip0;
#            aj = j*skip1;
#            al = l*skip2;
#
#
#            fprintf( output, "%e %e %e \n"
#                    ,g->k[0][i*skip0]
#                    ,g->k[1][j*skip1]
#                    ,zbcurva_c[ g->index( &ai, &aj, &al ) ]
#                    );
#
#
#        }
#
#    }
#
#    fflush(output);
#
#}



#void solidstructure::group_vel_vc_output( FILE *output, int skip0=1, int skip1=1, int skip2=1 )
#{
#
#
#
#    if (flag==0)
#    {
#
#        cout << "\n\n\n******....OJO...******\nExit, lattice_structure rutine need to be called!" << endl;
#        exit (EXIT_FAILURE);
#
#    }
#
#    for ( l = 0; l < g->N[2]/skip2; l++ )
#    for ( j = 0; j < g->N[1]/skip1; j++ )
#    {
#
#        fprintf( output, "\n");
#
#        for ( i = 0; i < g->N[0]/skip0; i++ )
#        {
#
#
#            //groupvel_c
#            ai = i*skip0;
#            aj = j*skip1;
#            al = l*skip2;
#
#
#            fprintf( output, "%e %e %e %e %e %e\n"
#                    ,g->k[0][i*skip0]
#                    ,g->k[1][j*skip1]
#                    ,groupvel_v[0][  g->index(&ai, &aj, &al) ]
#                    ,groupvel_v[1][  g->index(&ai, &aj, &al) ]
#                    ,groupvel_c[0][  g->index(&ai, &aj, &al) ]
#                    ,groupvel_c[1][  g->index(&ai, &aj, &al) ]
#                    );
#
#
#
#        }
#
#    }
#
#    fflush(output);
#}


#endif /* solidstructure_h */

