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

axesSym = ['x', 'y', 'z']
MaxDim = 3

# add options to see specific cross section of 3D data later
parser = argparse.ArgumentParser(description='Plot momentum matrix and energy dispersion.')
parser.add_argument('datapath', nargs='?', default='', type=str,help='path to folder contains data')
parser.add_argument('--noshow', dest='isShow', action='store_const', const=False, default=True, help = 'Do not show the plot (default: show)')
parser.add_argument('--f', dest='forceOverwrite', action='store_const', const=True, default=False, help = 'Force overwrite (default: ask)')

args = parser.parse_args()

BasicPath = os.getcwd()

FullPath = BasicPath + args.datapath

FileNameInput = FullPath + '/input.txt'
FileNamePmat = []
for i in range(MaxDim):
    FileNamePmat.append(FullPath + '/pmat' + axesSym[i] + '.bin')

FileNamePmatW = []
for i in range(MaxDim):
    FileNamePmatW.append(FullPath + '/pmat' + axesSym[i] + '_w.bin')

FileNamePmat = []
for i in range(MaxDim):
    FileNamePmat.append(FullPath + '/dmat' + axesSym[i] + '.bin')

FileNamePmatW = []
for i in range(MaxDim):
    FileNamePmatW.append(FullPath + '/dmat' + axesSym[i] + '_w.bin')

FileNameE  = FullPath + '/edispersion.txt'


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
plength = os.path.getsize(FileNamePmat[0]) / sizeofcomplex # assume double complex

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

# load energy dispersion
datE = np.loadtxt(FileNameE)    # (Nx*Ny*Nz, Nband)

print(f'Length of edispersion = ({len(datE[:,0])}, {len(datE[0,:])})')
print(f'Length of pmatx = {plength}')

# check length of datas
if (plength != Ntotal*(Nband*(Nband+1))/2 or len(datE[:,0]) != Ntotal or len(datE[0,:]) != Nband):
    print('Size of file is not corresponding to input file information')
    os.exit()
for i in range(dim):
    if (plength != os.path.getsize(FileNamePmat[i])/sizeofcomplex):
        print('Size of pmatx is not equal to + pmat' + axesSym[i])
        os.exit()

if (dim == 1):
    print('Sorry! plot with 1d data is not implemented yet!')
    os.exit()
if (dim == 3):
    print('Sorry! plot with 3d data is not implemented yet!')
    os.exit()

# load momentum matrix elements
print( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2 )
pmat = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    pmat[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmat[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

pmat_w = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    pmat_w[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmatW[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

dmat = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    dmat[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmat[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

dmat_w = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    dmat_w[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmatW[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

# now only 2D case are implemented!
# from here, implementation is 2d specific

Xgrid = np.zeros((Npoints[0], Npoints[1]))
Ygrid = np.zeros((Npoints[0], Npoints[1]))

for i in range(Npoints[0]):
    for j in range(Npoints[1]):
        Xgrid[i, j] = StartK[0] + float(i)/float(Npoints[0]-1) * Axes[0, 0] + float(j)/float(Npoints[1]-1) * Axes[1, 0]
        Ygrid[i, j] = StartK[1] + float(i)/float(Npoints[0]-1) * Axes[0, 1] + float(j)/float(Npoints[1]-1) * Axes[1, 1]

def gen_plots(data, name):
    #plt.figure()
    plt.pcolormesh(Xgrid, Ygrid, data, cmap=plt.get_cmap('seismic'))
    plt.colorbar()
    plt.savefig(name)
    plt.clf()

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            gen_plots(np.real(pmat[i, :, :, 0, index]), str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_r.png')
            gen_plots(np.imag(pmat[i, :, :, 0, index]), str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_i.png')
            
            gen_plots(np.real(pmat_w[i, :, :, 0, index]), str(m)+'_' + str(n)+'_pmat_w' + axesSym[i] + '_r.png')
            gen_plots(np.imag(pmat_w[i, :, :, 0, index]), str(m)+'_' + str(n)+'_pmat_w' + axesSym[i] + '_i.png')

            gen_plots(np.real(dmat[i, :, :, 0, index]), str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_r.png')
            gen_plots(np.imag(dmat[i, :, :, 0, index]), str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_i.png')
            
            gen_plots(np.real(dmat_w[i, :, :, 0, index]), str(m)+'_' + str(n)+'_dmat_w' + axesSym[i] + '_r.png')
            gen_plots(np.imag(dmat_w[i, :, :, 0, index]), str(m)+'_' + str(n)+'_dmat_w' + axesSym[i] + '_i.png')

#plt.show()
'''for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(pmat[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_r.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(pmat[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_i.png')
            plt.colorbar()
            plt.clf()
            index += 1

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(pmat_w[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_r_w.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(pmat_w[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_i_w.png')
            plt.colorbar()
            plt.clf()
            index += 1

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(pmat[i, :, :, 0, index] - pmat_w[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_r_diff.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(pmat[i, :, :, 0, index] - pmat_w[i, :, :, 0, index]))
            plt.savefig(str(m)+'_' + str(n)+'_pmat' + axesSym[i] + '_i_diff.png')
            plt.colorbar()
            plt.clf()
            index += 1

dmat = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    dmat[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmat[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

dmat_w = np.zeros(( dim, Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2), dtype=np.complex128)
for i in range(dim):
    dmat_w[i, :, :, :, :] = np.reshape(np.fromfile(FileNamePmatW[i], dtype=np.complex128), ( Npoints[0], Npoints[1], Npoints[2], (Nband*(Nband+1))//2) )

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(dmat[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_r.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(dmat[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_i.png')
            plt.colorbar()
            plt.clf()
            index += 1

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(dmat_w[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_r_w.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(dmat_w[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_i_w.png')
            plt.colorbar()
            plt.clf()
            index += 1

for i in range(dim):
    index = 0
    for m in range(Nband):
        for n in range(Nband):
            if (m > n):
                continue
            #p1 = plt.figure(1)
            plt.pcolormesh(Xgrid, Ygrid, np.real(dmat[i, :, :, 0, index] - dmat_w[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_r_diff.png')
            plt.colorbar()
            plt.clf()
            #p2 = plt.figure(2)
            plt.pcolormesh(Xgrid, Ygrid, np.imag(dmat[i, :, :, 0, index] - dmat_w[i, :, :, 0, index]), cmap=plt.get_cmap('seismic'))
            plt.savefig(str(m)+'_' + str(n)+'_dmat' + axesSym[i] + '_i_diff.png')
            plt.colorbar()
            plt.clf()
            index += 1'''