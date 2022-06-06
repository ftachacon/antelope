#!/usr/bin/env python3
"""
Reader and visual. of harmonic emission from solids, Chacon model """
import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm 
from matplotlib.ticker import MaxNLocator 
import numpy as np
import os 
import sys
import argparse
import shutil 
import cmath
import errno
from scipy.fftpack import fft, ifft

# add options to see specific cross section of 3D data later
parser = argparse.ArgumentParser(description='Plot momentum matrix and energy dispersion.')
parser.add_argument('datapath', nargs='?', default='', type=str,help='path to folder contains data')
parser.add_argument('--noshow', dest='isShow', action='store_const', const=False, default=True, help = 'Do not show the plot (default: show)')
parser.add_argument('--f', dest='forceOverwrite', action='store_const', const=True, default=False, help = 'Force overwrite (default: ask)')

args = parser.parse_args()

BasicPath = os.getcwd()

FullPath = BasicPath + "/" + args.datapath

FileNameInput = FullPath + '/input.txt'
FileNamePx = FullPath + '/pmatx.bin'
FileNamePy = FullPath + '/pmaty.bin'
FileNamePz = FullPath + '/pmatz.bin'
FileNameE  = FullPath + '/edispersion.txt'

axesSym = {'x', 'y', 'z'}
MaxDim = 3

Npoints = np.zeros(MaxDim, dtype=int)
for i in range(MaxDim):
    Npoints[i] = 1
Axes = np.zeros((MaxDim, MaxDim))
for i in range(MaxDim):
    Axes[i, i] = 1.0
StartK = np.zeros(MaxDim)

dim = -1
isStartKGiven = False
with open(FileNameInput) as f:
    lines = f.readlines()
    iline = 0
    while iline < len(lines):
        line1 = lines[iline].lstrip()
        if line1 and not line1.startswith('#'): # not empty and not starts with #
            line2 = line1.split()
            if (line2[0] == 'Npoints'):
                dim = len(line2) - 1
                for i in range(dim):
                    Npoints[i] = int(line2[i+1])
            elif (line2[0] == 'Nband'):
                Nband = int(line2[1])
            elif (line2[0] == 'Nvb'):
                Nvb = int(line2[1])
            elif (line2[0] == 'StartK'):
                if (dim < 0):
                    print('check input file order, Npoints must be faster than StartK')
                    os.exit()
                for i in range(dim):
                    StartK[i] = float(line2[i+1])
                isStartKGiven = True
            elif (line2[0] == 'Vectors'):
                if (dim < 0):
                    print('check input file order, Npoints must be faster than Vectors')
                    os.exit()
                itemp = 0
                while itemp < dim:
                    iline += 1
                    line1 = lines[iline].lstrip()
                    if line1 and not line1.startswith('#'):
                        line2 = line1.split()
                        for i in range(dim):
                            Axes[itemp][i] = float(line2[i])
                        itemp += 1
                if (not isStartKGiven):
                    for i in range(dim):
                        for j in range(dim):
                            StartK[i] -= Axes[j, i]/2.
        iline += 1

Ntotal = 1
for i in range(MaxDim):
    Ntotal = Ntotal * Npoints[i]
sizeofcomplex = 16
plength = os.path.getsize(FileNamePx) / sizeofcomplex # assume double complex
datE = np.loadtxt(FileNameE)

# print informations
print('Number of points = (', end='')
for i in range(dim):
    print(f' {Npoints[i]:4d}', end='')
print(' )')
print(f'Number of total points = {Ntotal}')
print(f'Number of bands = {Nband}')
print(f'Axes vectors:')
for i in range(dim):
    for j in range(dim):
        print(f'    {Axes[i, j]:.8f}', end='')
    print(' ')
print(f'StartK points:', end='')
for i in range(dim):
    print(f'    {StartK[i]:.8f}', end='')
print('')
print(f'Length of edispersion = ({len(datE[:,0])}, {len(datE[0,:])})')
print(f'Length of pmatx = {plength}')
if (plength != Ntotal*(Nband*(Nband+1))/2 or len(datE[:,0]) != Ntotal or len(datE[0,:]) != Nband):
    print('Size of file is not corresponding to input file information')
    os.exit()
if (dim > 1):
    if (plength != os.path.getsize(FileNamePy)/sizeofcomplex):
        print('Size of pmatx is not equal to pmaty')
        os.exit()
    if (dim > 2):
        if (plength != os.path.getsize(FileNamePz)/sizeofcomplex):
            print('Size of pmatx is not equal to pmatz')
            os.exit()
if (dim == 1):
    print('Sorry! plot with 1d data is not implemented yet!')
    os.exit()
if (dim == 3):
    print('Sorry! plot with 3d data is not implemented yet!')
    os.exit()


# now only 2D case are implemented!

'''pmat = np.zeros((dim, Npoints[0], Npoints[1], Npoints[2]))
Xgrid = np.zeros((Npoints[0], Npoints[1]))
Ygrid = np.zeros((Npoints[0], Npoints[1]))

for i in range(Npoints[0]):
    for j in range(Npoints[1]):
        Xgrid[i, j] = StartK[0] + float(i)/float(Npoints[0]-1) * Axes[0, 0] + float(j)/float(Npoints[1]-1) * Axes[1, 0]
        Ygrid[i, j] = StartK[1] + float(i)/float(Npoints[0]-1) * Axes[0, 1] + float(j)/float(Npoints[1]-1) * Axes[1, 1]

plt.pcolormesh(Xgrid, Ygrid, np.real())'''