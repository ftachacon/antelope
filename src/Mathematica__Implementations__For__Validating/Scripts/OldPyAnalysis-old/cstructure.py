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
Nvects  = 3;
Ndim    = 2;
Nsnap   = 160*2; #snaptshot per opt. cycle


###############################
#Lattice constant
a0       = 1./0.529;
t1       = 0.075;
t2       = t1/3.;
M0       = 2.54*t2;
phi0     = 1.16;
eps_phi0 = 0.5e-5;



###############################
## MOMENTUM AXIS
Nx      = 251; #51;   #250
Ny      = 100;   #100;  #101


xkmax   = 1.918999727359802;
ykmax   = 1.107935009166000;


dkx     = 2.*xkmax/( Nx - 1. );
dky     = 2.*ykmax/( Ny - 1. );



################################
## TIME STEPS
dt      = 1.401970981235;   #1.401970981234623283439
Ntime   = 5121;             #2560;             #
skiper  = int(Ntime/Nsnap);


w0      = 0.014;
T0      = 2.*pi/w0;
E0      = 0.007;


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
K[0]        = 4.*pi/np.sqrt(3.)/3./a0 ;
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
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])

    return plt


    
###########################
# Dipole matrix element
def inter_dipole(kx,ky,fx):
    integ = np.zeros(2);

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
    Nx  = leng(kx);
    dkx = kx[1]-kx[0];
    jv  = np.zeros(2);
    
    for n in range(Nx):
        gvelv = xygroup_velV2( kx[n], ky )
        gvelc = xygroup_velC2( kx[n], ky )
        
        jv[0]+= nc[n]*gvelc[0] + (1.-nc[n])*gvelv[0]
        jv[1]+= nc[n]*gvelc[1] + (1.-nc[n])*gvelv[1]

    jv*=    -dkx;
    return jv;





#####################################s
#####################################
#### anomalous velocity average  ####
def jva(kx,ky,nc,efield):
    # jvg[0] = nc*va_c[0] + nv*va_v[0]
    #
    Nx  = leng(kx);
    dkx = kx[1]-kx[0];
    ja  = np.zeros(2);

    for n in range(Nx):
        
        zcurvac0    =  zBerryCurvaC2( kx[n], ky  );
        zcurvav0    = -zcurvac0;
        janalc0     = -xyanomalous_velC( efield, zcurvac0 );
        janalv0     = -xyanomalous_velV( efield, zcurvav0 );
        
        ja[0]+= nc[n]*janalc0[0] + (1.-nc[n])*janalv0[0]
        ja[1]+= nc[n]*janalc0[1] + (1.-nc[n])*janalv0[1]
    
    ja*= -dkx;
    return ja;




#####################################
#####################################
## Intra band current calculation  ##
def j_ra(kx,ky,nc,efield):
    jgroupv0 = jvg(kx,ky,nc);
    janomal0 = jva(kx,ky,nc,efield);

    return jgroupv0+janomal0;




############################################
# Components of the Haldane model Hamiltonian
def B0( kx, ky ):

    b0res0 = 0.;
    
    
    for n in range(Nvects):
        b0res0+= np.cos( kx*vec_b[0,n] + ky*vec_b[1,n] );

    b0res0*= 2. * t2 * np.cos( phi0 );
    
    return b0res0;




######################################
##
def B1( kx, ky ):
    b1res0 = 0.;

    for n in range( Nvects ):
        b1res0+= np.cos( kx*vec_a[0,n] + ky*vec_a[1,n] );

    b1res0*= t1;

    return b1res0;


def B2( kx, ky ):
    b2res0 = 0.;

    for n in range( Nvects ):
        b2res0+= np.sin( kx*vec_a[0,n] + ky*vec_a[1,n] );

    b2res0*= t1;

    return b2res0;



def B3( kx, ky ):
    b3res0 = 0.;

    for n in range( Nvects ):
        b3res0+= np.sin( kx*vec_b[0,n] + ky*vec_b[1,n] );
    

    b3res0*= -2.0 * t2 * np.sin( phi0 );
    b3res0+= M0;

    return b3res0;

#
#
#
def BNorm( b1res0, b2res0, b3res0 ):


    bnorm0 = np.sqrt(
                   b1res0 * b1res0
                  +b2res0 * b2res0
                  +b3res0 * b3res0
                  );

    return bnorm0
#
#
#
def xgradB0( kx, ky ):


    xb0der0 = 0.;

    for n in range( Nvects ):
        xb0der0+=  vec_b[0,n]*np.sin( vec_b[0,n]*kx + vec_b[1,n]*ky );

    xb0der0*= -2.*t2*np.cos( phi0 );

    return xb0der0;



#
# GRADIENT ON THE HM HAMILTONIAN COMPONENTS
#

def ygradB0( kx, ky ):
    yb0der0 = 0.;

    for n in range( Nvects ):
        yb0der0+=  vec_b[1,n]*np.sin( vec_b[0,n]*kx + vec_b[1,n]*ky );

    yb0der0*= -2.*t2*np.cos( phi0 );

    return yb0der0


#

#
def xgradB1( kx, ky ):
    xb1der0 = 0.;

    for n in range(0, Nvects ):
        xb1der0+=  vec_a[0,n]*np.sin( vec_a[0,n]*kx + vec_a[1,n]*ky );

    xb1der0*= -t1;
    
    return xb1der0;


#
#
#
def ygradB1( kx, ky ):

    res = 0;

    for n in range( Nvects ):
        res+=  vec_a[1,n]*np.sin( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );

    res*= -t1;
    return res;

#
#
#
def xgradB2( kx, ky ):

    xb2der0 = 0.;

    for n in range( Nvects ):
        xb2der0+=  vec_a[0,n]*np.cos( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );


    xb2der0*= t1;
    return xb2der0

#
#
def ygradB2( kx, ky ):
    (res) = 0;

    for n in range( Nvects ):
        (res)+=  vec_a[1,n]*np.cos( vec_a[0,n]*(kx) + vec_a[1,n]*(ky) );


    (res)*= t1;
    return res
#
#
#
def xgradB3( kx, ky ):

    xb3der0 = 0;

    for n in range( Nvects ):
        xb3der0+=  vec_b[0,n]*np.cos( vec_b[0,n]*(kx) + vec_b[1,n]*(ky) );


    xb3der0*= -2.*t2*np.sin( phi0 );
    return xb3der0


#
#
#
def ygradB3( kx, ky ):
    (res) = 0.;

    for n in range( Nvects ):
        (res)+=  vec_b[1,n]*np.cos( vec_b[0,n]*(kx) + vec_b[1,n]*(ky) );

    (res)*= -2.*t2*np.sin( phi0 );
    return res;




def xgradBNorm( b1res0,  b2res0 ,b3res0
               ,xb1der0, xb2der0, xb3der0
               ,bnorm0 ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #
    #    xgradB1( &kx, &ky, &xb1der0 );
    #    xgradB2( &kx, &ky, &xb2der0 );
    #    xgradB3( &kx, &ky, &xb3der0 );
    #     */
    
    
    return (xb1der0*b1res0 + xb2der0*b2res0 + xb3der0*b3res0 )/bnorm0;



#
#
#
#
def ygradBNorm( b1res0, b2res0, b3res0, yb1der0, yb2der0, yb3der0, bnorm0  ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #
    #    ygradB1( &kx, &ky, &yb1der0 );
    #    ygradB2( &kx, &ky, &yb2der0 );
    #    ygradB3( &kx, &ky, &yb3der0 );
    #     */
    
    
    return ( yb1der0*b1res0 + yb2der0*b2res0 + yb3der0*b3res0 )/bnorm0;






def phiBloch(  b1res0, b2res0 ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #     B1( &kx, &ky );
    #     B2( &kx, &ky );
    #     */
    
    
    return np.arctan2( b2res0, b1res0 );



#
#
#
#
def thetaBloch( b3res0, bnorm0 ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #    B3( &kx, &ky );
    #    BNorm( &kx, &ky );
    #     */
    
    return np.arccos( b3res0 / bnorm0 );


#
#
#
#
def xgrad_rphiBloch( b1res0, b2res0,  xb1der0, xb2der0  ):
    #    /*
    #     NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #
    #     B1( &kx, &ky );
    #     B2( &kx, &ky );
    #
    #     xgradB1( &kx, &ky, &xb1der0 );
    #     xgradB2( &kx, &ky, &xb2der0 );
    #     */
    
    return   ( b1res0 * xb2der0 - b2res0 * xb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0 + eps_phi0 );


#
#
#
#
def xgrad_phiBloch( b1res0, b2res0, xb1der0, xb2der0 ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #
    #     B1( &kx, &ky );
    #     B2( &kx, &ky );
    #
    #     xgradB1( &kx, &ky, &xb1der0 );
    #     xgradB2( &kx, &ky, &xb2der0 );
    #
    #    */
    
    return ( b1res0 * xb2der0 - b2res0 * xb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0 );






def ygrad_rphiBloch( b1res0, b2res0, yb1der0, yb2der0 ):
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate,
    #
    #    B1( &kx, &ky );
    #    B2( &kx, &ky );
    #
    #    ygradB1( &kx, &ky, &yb1der0 );
    #    ygradB2( &kx, &ky, &yb2der0 );
    #
    #     */
    
    return ( b1res0 * yb2der0 - b2res0 * yb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0 + eps_phi0 );


#
#
#
#
def ygrad_phiBloch( b1res0,b2res0,yb1der0, yb2der0  ):
    
    
    #    /* NOTE
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #     B1( &kx, &ky );
    #     B2( &kx, &ky );
    #
    #     ygradB1( &kx, &ky, &yb1der0 );
    #     ygradB2( &kx, &ky, &yb2der0 );
    #
    #     */
    
    return ( b1res0 * yb2der0 - b2res0 * yb1der0 )/( b1res0 * b1res0 + b2res0 * b2res0  );


#
#
#
def xgrad_thetaBloch(  b3res0, bnorm0, xb3der0, xbnorm0  ):
    #{
    
    #    /*
    #     Before using this fungtion we should called and evaluate:
    #
    #     B3( &kx, &ky );
    #     BNorm( &kx, &ky );
    #
    #
    #     xgradB3( &kx, &ky, &xb3der0 );
    #     xgradBNorm( &kx, &ky, &xbnorm0 );
    #
    #    */
    temp =  bnorm0 * bnorm0 - b3res0 * b3res0
    #print( "\nTemDenoValue0 = ",temp)
    if temp == 0.:
        print("***************\nOJO... denominator is zero = ",temp)
        print("But, it is regularized by eps = ",eps_phi0,"\n*******\n")
    #temp = eps_phi0;
        #sys.exit()
    
    return -( bnorm0 * xb3der0 -  xbnorm0 * b3res0 ) / np.sqrt( temp )/ bnorm0 ;




#}
#
#
def xgrad_rthetaBloch( b3res0, bnorm0, xb3der0, xbnorm0 ):
    #{
    #
    #    /*
    #     Before using this fungtion we should called and evaluate:
    #
    #     B3( &kx, &ky );
    #     BNorm( &kx, &ky );
    #
    #
    #     xgradB3( &kx, &ky, &xb3der0 );
    #     xgradBNorm( &kx, &ky, &xbnorm0 );
    #
    #     */
    
    
    return ( bnorm0 * xb3der0 -  xbnorm0 * b3res0 )/( np.sqrt( bnorm0 * bnorm0 - b3res0 * b3res0  ) + eps_phi0 )/( bnorm0 + eps_phi0 );

#}
#
#
def ygrad_thetaBloch( b3res0, bnorm0, yb3der0, ybnorm0 ):
    #{
    #
    #
    #    /*
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #     B3( &kx, &ky );
    #     BNorm( &kx, &ky );
    #
    #
    #     xgradB3( &kx, &ky, &yb3der0 );
    #     xgradBNorm( &kx, &ky, &ybnorm0 );
    #
    #     */
    
    temp =  bnorm0 * bnorm0 - b3res0 * b3res0
    #print( "\nTemDenoValue0 = ",temp)
    if temp == 0.:
        print("***************\nOJO... denominator is zero = ",temp)
        print("But, it is regularized by eps = ",eps_phi0,"\n*******\n")
        temp =eps_phi0;
    
    
    return -( bnorm0 * yb3der0 -  ybnorm0 * b3res0 )/ np.sqrt( temp  ) / bnorm0  ;






def ygrad_rthetaBloch( b3res0, bnorm0, yb3der0, ybnorm0 ):
    #{
    
    #
    #    /*
    #
    #     Before using this fungtion we should called and evaluate:
    #
    #     B3( &kx, &ky );
    #     BNorm( &kx, &ky );
    #
    #
    #     xgradB3( &kx, &ky, &yb3der0 );
    #     xgradBNorm( &kx, &ky, &ybnorm0 );
    #
    #    */
    
    return -( bnorm0 * yb3der0 -  ybnorm0 * b3res0 )/( np.sqrt( bnorm0 * bnorm0 - b3res0 * b3res0  ) + eps_phi0 )/( bnorm0 + eps_phi0 );






####################################################
####################################################
#
# ENERGY DISPERSIONS CONDUCTION AND VALENCE BANDS
def Ec( b0res0, bnorm0 ):
    return b0res0 + bnorm0;


def Ec2(kx,ky):
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients(kx,ky);
    return b0res0 + bnorm0;


#
#
def Ev( b0res0, bnorm0 ):
    return b0res0 - bnorm0;


def Ev2(kx,ky):
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients(kx,ky);
    return b0res0 - bnorm0;

#
#
#
def Eg( bnorm0 ):
    return 2.*bnorm0;


def Eg2( kx, ky ):
    return Ec2(kx,ky)-Ev2(kx,ky)







####################################################
#####################################################
## Conduction Band Berry Connection
def xBerryConnectionC( b3res0, bnorm0, xphi_rbloch0 ):
    return b3res0/2./( bnorm0 + eps_phi0)*xphi_rbloch0;





def xBerryConnectionC2( kx, ky  ):
    b1res0 = B1( kx, ky );
    b2res0 = B2( kx, ky );
    b3res0 = B3( kx, ky );
    bnorm0 = BNorm( b1res0, b2res0, b3res0 );
    
    xb1der0 = xgradB1( kx, ky );
    xb2der0 = xgradB2( kx, ky );
    xphi_rbloch0     = xgrad_rphiBloch( b1res0, b2res0, xb1der0, xb2der0 );
    
    return b3res0/2./( bnorm0 + eps_phi0)*xphi_rbloch0;



def xBerryConnectionV( b3res0, bnorm0, xphi_rbloch0 ):
    #    //This function requires the evaluation of B3(...), BNorm(...) and xgrad_phiBloch before being called
    return -b3res0/2./( bnorm0 + eps_phi0)*xphi_rbloch0;



def xBerryConnectionV2( kx,ky  ):
    
    b1res0 = B1( kx, ky );
    b2res0 = B2( kx, ky );
    b3res0 = B3( kx, ky );
    bnorm0 = BNorm( b1res0, b2res0, b3res0 );
    
    xb1der0 = xgradB1( kx, ky );
    xb2der0 = xgradB2( kx, ky );
    xphi_rbloch0     = xgrad_rphiBloch( b1res0, b2res0, xb1der0, xb2der0 );
    
    return -b3res0/2./( bnorm0 + eps_phi0)*xphi_rbloch0;



def yBerryConnectionC( b3res0, bnorm0, yphi_rbloch0 ):
    return b3res0/2./( bnorm0 + eps_phi0)*yphi_rbloch0;




def yBerryConnectionC2( kx,ky ):
    b1res0 = B1( kx, ky );
    b2res0 = B2( kx, ky );
    b3res0 = B3( kx, ky );
    bnorm0 = BNorm( b1res0, b2res0, b3res0 );
    
    yb1der0 = ygradB1( kx, ky );
    yb2der0 = ygradB2( kx, ky );
    yphi_rbloch0     = ygrad_rphiBloch( b1res0, b2res0, yb1der0, yb2der0 );
    
    return b3res0/2./( bnorm0 + eps_phi0)*yphi_rbloch0;






def yBerryConnectionV( b3res0, bnorm0, yphi_rbloch0 ):
#    //This function requires the evaluation of B3(...), BNorm(...) and ygrad_phiBloch before being called
    return -b3res0/2./( bnorm0 + eps_phi0)*yphi_rbloch0;





def yBerryConnectionV2( kx, ky ):
    b1res0 = B1( kx, ky );
    b2res0 = B2( kx, ky );
    b3res0 = B3( kx, ky );
    bnorm0 = BNorm( b1res0, b2res0, b3res0 );
    
    yb1der0 = ygradB1( kx, ky );
    yb2der0 = ygradB2( kx, ky );
    yphi_rbloch0     = ygrad_rphiBloch( b1res0, b2res0, yb1der0, yb2der0 );
    
    return -b3res0/2./( bnorm0 + eps_phi0)*yphi_rbloch0;







#//Berry Connection difference
def xyChig( b3res0, bnorm0, xphi_rbloch0, yphi_rbloch0 ):

    xcberryconn0 = 2.0 * xBerryConnectionC( b3res0, bnorm0, xphi_rbloch0 );
    ycberryconn0 = 2.0 * yBerryConnectionC( b3res0, bnorm0, yphi_rbloch0 );


    return [xcberryconn0,ycberryconn0];




def xyChig2( kx,ky ):
    xcberryconn0 = 2.0 * xBerryConnectionC2( kx,ky );
    ycberryconn0 = 2.0 * yBerryConnectionC2( kx,ky );
    
    return [xcberryconn0,ycberryconn0];




####################################################
####################################################
#//Dipole Matrix Element from Val. to Cond. Band
def xydipoleCV( theta_bloch0
             ,xphi_bloch0
             ,yphi_bloch0
             ,xtheta_bloch0
             ,ytheta_bloch0):


    xdip_cv0 =  np.sin( theta_bloch0 ) * xphi_bloch0/2. + I * xtheta_bloch0/2.;
    ydip_cv0 =  np.sin( theta_bloch0 ) * yphi_bloch0/2. + I * ytheta_bloch0/2.;


    return [xdip_cv0,ydip_cv0]


##
def xydipoleCV2( kx, ky):
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients(kx,ky)
    xb0der0, xb1der0, xb2der0, xb3der0          = set_of_xbgrad( kx, ky );
    yb0der0,yb1der0,yb2der0,yb3der0             = set_of_ybgrad( kx, ky );
    
    xbnorm0         = xgradBNorm( b1res0, b2res0, b3res0, xb1der0 ,xb2der0, xb3der0, bnorm0 );
    ybnorm0         = ygradBNorm( b1res0, b2res0, b3res0, yb1der0 ,yb2der0, yb3der0, bnorm0 );
    
    theta_bloch0    = thetaBloch( b3res0, bnorm0 );
    
    xphi_bloch0     = xgrad_phiBloch( b1res0, b2res0, xb1der0, xb2der0 );
    yphi_bloch0     = ygrad_phiBloch( b1res0, b2res0, yb1der0, yb2der0 );
                                                  
    xtheta_bloch0   = xgrad_thetaBloch( b3res0, bnorm0, xb3der0, xbnorm0 );
    ytheta_bloch0   = ygrad_thetaBloch( b3res0, bnorm0, yb3der0, ybnorm0 );
    
    xdip_cv0        =  np.sin( theta_bloch0 ) * xphi_bloch0/2. + I * xtheta_bloch0/2.;
    ydip_cv0        =  np.sin( theta_bloch0 ) * yphi_bloch0/2. + I * ytheta_bloch0/2.;
    
    
    return [xdip_cv0,ydip_cv0]


####################################################
####################################################
# GROUP VELOCITIES
# conduction band group vel.
def xygroup_velC( xb0der0, yb0der0, xbnorm0, ybnorm0 ):
    xgroup_c0 =  xb0der0 + xbnorm0 ;
    ygroup_c0 =  yb0der0 + ybnorm0 ;

    return [xgroup_c0, ygroup_c0]



def set_of_bcohefficients(kx,ky):
    b0res0 = B0( kx, ky );
    b1res0 = B1( kx, ky );
    b2res0 = B2( kx, ky );
    b3res0 = B3( kx, ky );
    bnorm0 = BNorm( b1res0, b2res0, b3res0 );

    return [b0res0,b1res0,b2res0,b3res0,bnorm0]



def set_of_xbgrad(kx,ky):
    xb0der0 = xgradB0( kx, ky );
    xb1der0 = xgradB1( kx, ky );
    xb2der0 = xgradB2( kx, ky );
    xb3der0 = xgradB3( kx, ky );
    return [xb0der0,xb1der0,xb2der0,xb3der0]



def set_of_ybgrad(kx,ky):
    yb0der0 = ygradB0( kx, ky );
    yb1der0 = ygradB1( kx, ky );
    yb2der0 = ygradB2( kx, ky );
    yb3der0 = ygradB3( kx, ky );
    return [yb0der0,yb1der0,yb2der0,yb3der0]



def xygroup_velC2( kx, ky ):
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients(kx,ky)
    xb0der0, xb1der0, xb2der0, xb3der0          = set_of_xbgrad( kx, ky );
    yb0der0,yb1der0,yb2der0,yb3der0             = set_of_ybgrad( kx, ky );
    
    xbnorm0 = xgradBNorm( b1res0, b2res0, b3res0, xb1der0 ,xb2der0, xb3der0, bnorm0 );
    ybnorm0 = ygradBNorm( b1res0, b2res0, b3res0, yb1der0 ,yb2der0, yb3der0, bnorm0 );

    xgroup_c0 =  xb0der0 + xbnorm0 ;
    ygroup_c0 =  yb0der0 + ybnorm0 ;
    return [xgroup_c0, ygroup_c0]



def xygroup_velV( xb0der0, yb0der0, xbnorm0, ybnorm0 ):
    xgroup_v0 =  xb0der0 - xbnorm0 ;
    ygroup_v0 =  yb0der0 - ybnorm0 ;
    return [xgroup_v0,ygroup_v0]



def xygroup_velV2( kx, ky ):
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients(kx,ky)
    xb0der0, xb1der0, xb2der0, xb3der0          = set_of_xbgrad( kx, ky );
    yb0der0,yb1der0,yb2der0,yb3der0             = set_of_ybgrad( kx, ky );
    
    
    xbnorm0 = xgradBNorm( b1res0, b2res0, b3res0, xb1der0 ,xb2der0, xb3der0, bnorm0 );
    ybnorm0 = ygradBNorm( b1res0, b2res0, b3res0, yb1der0 ,yb2der0, yb3der0, bnorm0 );
    
    
    #print ("\nybnorm0 = ",ybnorm0)
    
    xgroup_v0 =  xb0der0 - xbnorm0 ;
    ygroup_v0 =  yb0der0 - ybnorm0 ;
    
    return [xgroup_v0,ygroup_v0]








####################################################
####################################################
##  BERRY CURVATUREs
#//Conduction band Berry Curvature: -((xGradtheta*yGradphi - xGradphi*yGradtheta)*Sin(theta))/2.
def zBerryCurvaC( theta_bloch0
                 ,xphi_bloch0
                 ,yphi_bloch0
                 ,xtheta_bloch0
                 ,ytheta_bloch0 ):
    #/*
    # Before using this fungtion we should called:
    #
    #    thetaBloch( &kx, &ky, &theta_bloch0 );
    #
    #
    #        xgrad_thetaBloch( &kx, &ky, &xtheta_bloch0 );
    #        xgrad_phiBloch( &kx, &ky, &xphi_bloch0 );
    #
    #
    #        ygrad_thetaBloch( &kx, &ky, &ytheta_bloch0 );
    #        ygrad_phiBloch( &kx, &ky, &yphi_bloch0 );
    #
    # */
    
    
    
    return - 1./2. * np.sin( theta_bloch0  ) * ( xtheta_bloch0 * yphi_bloch0 - xphi_bloch0 * ytheta_bloch0 );



def zBerryCurvaC2( kx, ky  ):
    
    b0res0,  b1res0,  b2res0, b3res0, bnorm0    = set_of_bcohefficients( kx, ky );
    xb0der0, xb1der0, xb2der0, xb3der0          = set_of_xbgrad( kx, ky );
    yb0der0,yb1der0,yb2der0,yb3der0             = set_of_ybgrad( kx, ky );
    
    theta_bloch0 = thetaBloch( b3res0, bnorm0 );
    
    xbnorm0 = xgradBNorm( b1res0, b2res0, b3res0, xb1der0, xb2der0, xb3der0, bnorm0 );
    ybnorm0 = ygradBNorm( b1res0, b2res0, b3res0, yb1der0, yb2der0, yb3der0, bnorm0 );
    
    xphi_bloch0     = xgrad_phiBloch( b1res0, b2res0, xb1der0, xb2der0 );
    yphi_bloch0     = ygrad_phiBloch( b1res0, b2res0, yb1der0, yb2der0 );
    
    xtheta_bloch0   = xgrad_thetaBloch( b3res0, bnorm0, xb3der0, xbnorm0 );
    ytheta_bloch0   = ygrad_thetaBloch( b3res0, bnorm0, yb3der0, ybnorm0 );
    
    return - 1./2. * np.sin( theta_bloch0  ) * ( xtheta_bloch0 * yphi_bloch0 - xphi_bloch0 * ytheta_bloch0 );



def zBerryCurvaV( theta_bloch0
                 ,xphi_bloch0
                 ,yphi_bloch0
                 ,xtheta_bloch0
                 ,ytheta_bloch0 ):
    
    
    
    return -zBerryCurvaC( theta_bloch0, xphi_bloch0, yphi_bloch0, xtheta_bloch0, ytheta_bloch0 )







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




####################################################
####################################################
# Chern No. calculation
def Chern_Number(kx,ky):
    Nx  = len(kx);
    Ny  = len(ky);
    dkx = kx[1]-kx[0];
    dky = ky[1]-ky[0];

    chernNo = 0.;
    for j in range(Ny):
        for i in range(Nx):
            zbcuva = zBerryCurvaC2( kx[i], ky[j]  );
            chernNo+= zbcuva;

    return chernNo*dkx*dky/4./pi

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
