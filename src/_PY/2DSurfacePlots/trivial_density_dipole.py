#!/usr/bin/python
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import numpy as np
import os
import sys
import shutil



import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker

from matplotlib.pylab import *
from numpy import arange
from scipy.interpolate import interp2d

set_of_params   = sys.argv;
iparam             =  set_of_params[1];  #mag. flux phase

#################################
#
kx 	 = np.loadtxt('kxAxis.dat');
ky   = np.loadtxt('kyAxis.dat');

temp = 'RealDipoleTrivialGauge'+iparam+'.dat';

Dip1 = np.transpose( np.loadtxt( temp ) );
print '\nSize of Dip1 = ', np.shape( Dip1 );

X,Y     = np.meshgrid(kx,ky);

nmin    = -2#Dip1.min();
nmax    = 2.01#Dip1.max();

#################################
width   = 11;
hight   = width#/1.62;
fig     = plt.figure( figsize=(width,hight) );
ax0     = plt.axes([0.11, 0.115, 0.865, 0.87]);
for axis1 in ['top','bottom','left','right']:
    ax0.spines[axis1].set_linewidth(2.5);

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );
colorMap0='seismic'#'CMRmap'#'RdBu_r'#'gnuplot2'#'bwr'#'coolwarm'#'gist_rainbow'#
cmap = plt.get_cmap(colorMap0);
norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );

f = interp2d(kx, ky, Dip1, kind='cubic');
xnew = np.arange( min(kx), max(kx), max(kx)*2/700 )
ynew = np.arange( min(ky), max(ky), max(ky)*2/700)
Xn, Yn = np.meshgrid(xnew, ynew)
Dat2    = f(xnew,ynew)

cf = plt.pcolormesh( Xn, Yn, Dat2,
                 vmin=nmin,
                 vmax=nmax,
                 cmap = cmap );

plt.tick_params(labelsize=23);
xticks0  = np.arange(-2,2,1);
plt.xticks( xticks0 );
plt.xlim([min(kx),max(kx)]);
plt.ylim([min(ky),max(ky)]);

K1 = np.array([-1.279333,0])
K2 = np.array([1.279333,0])
K3 = np.array([-1.279333/2,1.1079])
K4 = np.array([1.279333/2,1.1079])
K5 = np.array([-1.279333/2,-1.1079])
K6 = np.array([1.279333/2,-1.1079])

plt.plot([K1[0],K3[0]],[K1[1],K3[1]],'g',lw=3)
plt.plot([K4[0],K2[0]],[K4[1],K2[1]],'g',lw=3)
plt.plot([K6[0],K2[0]],[K6[1],K2[1]],'g',lw=3)
plt.plot([K3[0],K4[0]],[K3[1],K4[1]],'g',lw=3)
plt.plot([K1[0],K5[0]],[K1[1],K5[1]],'g',lw=3)
plt.plot([K5[0],K6[0]],[K5[1],K6[1]],'g',lw=3)

wid=21
plt.plot(K1[0],K1[1],'wo',lw=6,markersize=wid)
plt.plot(K6[0],K6[1],'wo',lw=6,markersize=wid)
plt.plot(K4[0],K4[1],'wo',lw=6,markersize=wid)

plt.plot(K5[0],K5[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K3[0],K3[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)
plt.plot(K2[0],K2[1],color=[1,0.9,0],marker='o',lw=6,markersize=wid)

plt.xlabel(r'$\rm k_x\, (a.u.)$', fontsize=33);
if iparam=='1':
    plt.ylabel(r'$\rm k_y\, (a.u.)$', fontsize=33);
    plt.text(-1.82, 2.0, r'$\rm (a) \,\, Gauge\, A$', fontsize=30, color='k')
    plt.text(-1.82, 1.50, r"$\rm  Trivial\,\, Dipole $", fontsize=26, color='k')

if iparam=='2':
    #plt.ylabel(r'$\rm k_y\, (a.u.)$', fontsize=33);
    plt.text(-1.82, 2.0, r'$\rm (b) \,\, Gauge\, B$', fontsize=30, color='k')
    plt.text(-1.82, 1.50, r"$\rm  Trivial\,\, Dipole $", fontsize=26, color='k')

fname               = 'DipoleGauge'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig( fileNamePicture, dpi = 150 );

plt.show()
