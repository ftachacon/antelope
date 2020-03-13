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
n 	     = 0;
ntemp    = 34;
Nt 	     = 133411;
Nthalf 	 = int(133411/16./2.);

if iparam=="LCP":
    f 	     = open("files1.dat","r");
else:
    f        = open("files2.dat","r");

myinfo   = [];
phase 	 = [];
data 	 = [];
hhgMAP0  = [];
hhgMAP   = np.zeros( (ntemp,Nthalf) )#[];


#################################
for line in f:
    myinfo.append(line);
    data.append(  np.loadtxt( line.strip() ) );
    phase.append( data[n][0,1] );
    hhgMAP0.append( data[n][1:Nthalf+1,5] );
    n+=1;
    print("\nfile = ",line," n = ",n);

hhgMAP   = np.zeros( (n,Nthalf) )#[];
for i in range(n):
    hhgMAP[i,:] = hhgMAP0[i][:]

om = data[0][1:Nthalf+1,0];

print( len(data[0][:,0]), 'omega-max =  ', om[-1] );
print( phase );
print( '\nNlist = ', n, '  Nthalf = ', Nthalf );

HHG_MAP = hhgMAP #np.reshape( hhgMAP, (n,Nt-1) );



#################################
width   = 11;
hight   = width/1.62;

fig     = plt.figure( figsize=(width,hight) );
#ax1     = fig.add_axes([0.10, 0.125, 0.92, 0.825])
plt.axes([0.07, 0.125, 0.79, 0.825])
plt.rc('lines', lw=3, color='k')

#nmin    = np.log10(1e-16);
#nmax    = np.log10(1e-2);

nmin    = 1e-16;
nmax    = 1e-1;

levels = MaxNLocator(nbins=28).tick_values( nmin, nmax );


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
colorMap0='CMRmap'#'RdBu_r'#'gnuplot2'#'bwr'#'coolwarm'#'gist_rainbow'#    #'seismic'#'gray_r'#'RdGy'#'PRGn'#'bwr'#'gnuplot2'#'CMRmap_r'#'coolwarm'#
cmap = plt.get_cmap(colorMap0);

norm = BoundaryNorm(levels,ncolors=cmap.N,clip=True );
X,Y     = np.meshgrid(om,phase);

print ('shape of X ',X.shape);

f       = interp2d(om, phase, np.log10(HHG_MAP+1e-17), kind='cubic')
xnew = np.arange( 0., max(om), float(om[1]-om[0]) )
ynew = np.arange( 0., max(phase), .01)

data1 =  f(xnew,ynew)
#data1=np.log10(HHG_MAP+1e-17)
Xn, Yn = np.meshgrid(xnew, ynew)


#pcolormesh contourf
#cf  = plt.pcolormesh( X, Y, np.power(10.,data1),
#                norm  = LogNorm( vmin=nmin, vmax=nmax ),
#                 cmap = cmap );

cf = plt.pcolormesh( Xn, Yn, np.power(10.,data1),
                norm  = LogNorm( vmin=nmin, vmax=nmax ),
                 cmap = cmap );
#cf = plt.pcolormesh( Xn, Yn, data1
#                    ,vmin=nmin
#                    ,vmax=nmax #,levels=levels
#                    ,cmap=cmap );


xticks0  = np.arange(1,32,4);
plt.xticks( xticks0 );
plt.xlabel(r'$\rm Harmonic\,\, Order$', fontsize=31);

if iparam=='LCP':
    cb = plt.colorbar( cf );#mappable
    cb.ax.tick_params(labelsize=25);

plt.tick_params(labelsize=25);
plt.ylim([0,3.12]);
plt.xlim([0, 33.]);


if iparam=='LCP':
    plt.ylabel(r'$\phi_0 {\rm (rad.)}$', fontsize=31);
else:
    plt.ylabel("", fontsize=31);

if iparam=='RCP':
    chernNo = np.loadtxt('Phi__Vs__Chern__Num.dat');
    gapPhi = np.loadtxt('Phi__Vs__Energy__Gap.dat');
    
    ax1=plt.axes([0.865, 0.125, 0.112, 0.825]);
    #plt.rc('lines', lw=3, color='k')
    ax1.plot(gapPhi[:,1],gapPhi[:,0],lw=2);
    ax1.set_yticks([]);
    ax1.set_xticks([0,0.1])#,('0','0.1'));
    ax1.set_xlim([0.0,max(gapPhi[:,1])]);
    ax1.set_ylim([0,3.12]);
    ax1.tick_params(labelsize=20);
    for axis in ['top','bottom','left','right']:
        ax1.spines[axis].set_linewidth(2.5)



#cb.set_tick_params(labelsize=16)
#plt.xlim([0, 27])


fname = 'HHG__AND__MAP__'+iparam+'.png'
filename1           =   fname;
fileNamePicture     =   filename1;
plt.savefig(fileNamePicture, dpi = 150);

#ntotal = ntemp*Nthalf
#outdata=np.zeros( (ntotal,3))
#outname='DATA__HHG__PHASE__MAP__LCP.dat'
outname='DATA__HHG__PHASE__MAP__' + iparam + '.dat'
afile   = open( outname, 'w' );

for i in range(n):
    for j in range(Nthalf):
        temp   = str('%.4e'%phase[i]) + '     ' + str( '%.4e'%om[j] ) +'     ' + str( '%.16e'%HHG_MAP[i,j] );
        afile.write( temp );
        afile.write("\n")
afile.close()


outname='DATA__PHASE__AXIS__' + iparam + '.dat'
afile   = open( outname, 'w' );
for i in range(n):
    temp   = str('%.4e'%phase[i]) ;
    afile.write( temp );
    afile.write("\n")
afile.close()

outname='DATA__FREQ__AXIS__' + iparam + '.dat'
afile   = open( outname, 'w' );
for j in range(Nthalf):
    temp   =  str( '%.4e'%om[j] )
    afile.write( temp );
    afile.write("\n")
afile.close()



plt.show()
