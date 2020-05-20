#! /usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

# input example
# ./polarizatoin_analysis.py dataname Ntheta order1 order2 ...
set_of_params   = sys.argv
argc = len(set_of_params)
dataname = set_of_params[1]
filename = dataname[:-4] 
#Ntheta = int(set_of_params[2])
#thetalist = np.linspace(0, np.pi, Ntheta)

#orderlist = np.zeros(argc-3, dtype=np.int32)
#Norder = argc-3
#for i in range(Norder):
#    orderlist[i] = int(set_of_params[i+3])
data = np.loadtxt(dataname)
elist = np.zeros(len(data[:,0])-1)
thetalist = np.zeros(len(data[:,0])-1)
rectphaselist = np.zeros(len(data[:,0])-1)

spectrum_xmin = 1
spectrum_xmax = 20
integrateWidth = 0.1
orderlist = np.array(range(spectrum_xmin, spectrum_xmax+1))
spectrum_xmin -= 1
spectrum_xmax += 1
roffset = 2.0
Norder = len(orderlist)
orderJlist = np.zeros((Norder, 2), dtype=np.complex128)
orderThetaList = np.zeros(Norder)
orderEList = np.zeros(Norder)
orderRectPhaseList = np.zeros(Norder)

xticks0  = np.arange(-4,20,1)

width = 11
hight = width/1.62
for j in range(1, len(data[:,0])):
    if (data[j, 0] < spectrum_xmin):
        continue
    if (data[j, 0] > spectrum_xmax):
        break
    rectphaselist[j-1] = 0.5 * np.angle((data[j, 1] + 1j*data[j,2])**2 + (data[j,3] + 1j*data[j,4])**2)
    newex = np.exp(-1j*rectphaselist[j-1]) * (data[j,1] + 1j*data[j,2])
    newey = np.exp(-1j*rectphaselist[j-1]) * (data[j,3] + 1j*data[j,4])
    a = np.abs(np.real(newex)+1j*np.real(newey))
    b = np.abs(np.imag(newex)+1j*np.imag(newey))
    #if (a > b):
    elist[j-1] = b/a
    if (np.real(newex) < 0):
        rectphaselist[j-1] += np.pi
        newex = -newex
        newey = -newey
    thetalist[j-1] = np.arctan2(np.real(newey)/a, np.real(newex)/a )

# integration for each order
dw = data[2, 0] - data[1, 0]
for j in range(1, len(data[:,0])):
    if (data[j, 0] < spectrum_xmin):
        continue
    if (data[j, 0] > spectrum_xmax):
        break
    for i in range(Norder):
        if (abs(data[j,0] - orderlist[i]) < integrateWidth):
            orderJlist[i, 0] += (data[j,1] + 1j*data[j,2]) * dw
            orderJlist[i, 1] += (data[j,3] + 1j*data[j,4]) * dw

for i in range(Norder):
    orderRectPhaseList[i] = 0.5 * np.angle((orderJlist[i,0])**2 + (orderJlist[i,1])**2)
    newex = np.exp(-1j*orderRectPhaseList[i]) * orderJlist[i, 0]
    newey = np.exp(-1j*orderRectPhaseList[i]) * orderJlist[i, 1]
    a = np.abs(np.real(newex)+1j*np.real(newey))
    b = np.abs(np.imag(newex)+1j*np.imag(newey))
    orderEList[i] = b/a
    if (np.real(newex) < 0):
        orderRectPhaseList[i] += np.pi
        newex = -newex
        newey = -newey
    orderThetaList[i] = np.arctan2(np.real(newey)/a, np.real(newex)/a )

fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
#plt.plot(data[1:,0], elist, '-bo')
plt.plot(data[1:,0], elist)
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel('ellipticity (abs)')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_e_raw.png')

fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
#plt.plot(data[1:,0], np.unwrap(thetalist)/np.pi, '-bo')
plt.plot(data[1:,0], np.unwrap(thetalist)/np.pi)
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel(r'$\theta/\pi$ (major-axis angle)')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_theta_raw.png')

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
ax.scatter(thetalist, data[1:,0], s=5.0)
ax.set_rorigin(-roffset)
#ax.set_rticks(orderlist)
ax.set_xticks(np.linspace(0, 2*np.pi-np.pi/4, 8))
plt.ylim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_theta_raw_polar.png')

fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
#plt.plot(data[1:,0], np.unwrap(rectphaselist)/np.pi, '-bo')
plt.plot(data[1:,0], np.unwrap(rectphaselist)/np.pi)
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel(r'$\gamma/\pi$ (rectifying angle)')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_gamma_raw.png')

###############
fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
plt.plot(orderlist, orderEList, '-bo')
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel('ellipticity (abs) - integrated')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_e.png')

fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
plt.plot(orderlist, np.unwrap(orderThetaList)/np.pi, '-bo')
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel(r'$\theta/\pi$ (major-axis angle) - integrated')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_theta.png')

fig  = plt.figure(figsize=(width,hight) )
ax1  = fig.add_axes([0.2, 0.15, 0.75, 0.75])
plt.plot(orderlist, np.unwrap(orderRectPhaseList)/np.pi, '-bo')
plt.grid(True)
plt.xlabel('Harmonic order')
plt.ylabel(r'$\gamma/\pi$ (rectifying angle) - integrated')
plt.xticks( xticks0 )
plt.xlim(spectrum_xmin, spectrum_xmax)
plt.savefig(filename + '_gamma.png')

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
#for i in range(Norder):
ax.scatter(orderThetaList, orderlist)
ax.set_rorigin(-roffset)
#ax.set_rticks(orderlist)
ax.set_xticks(np.linspace(0, 2*np.pi-np.pi/4, 8))
plt.savefig(filename + '_theta_polar.png')

plt.show()

# find index for each orders
'''xylist = np.zeros((Norder, 2), dtype=np.complex128)
for i in range(Norder):
    data = np.loadtxt(dataname)
    for j in range(1, len(data[:,0])):
        if (data[j, 0] > orderlist[i] ):
            if ((data[j-1, 0] + data[j, 0])/2 > orderlist[i]):
                index = j-1
            else:
                index = j
            break
    xylist[i, 0] = data[index, 1] + 1j*data[index, 2]
    xylist[i, 1] = data[index, 3] + 1j*data[index, 4]

# calculate index for each angles
outlist = np.zeros((Norder, Ntheta))
for i in range(Norder):
    for j in range(Ntheta):
        outlist[i, j] = abs(xylist[i, 0]*np.cos(thetalist[j]) + xylist[i, 1]*np.sin(thetalist[j]))**2

# normalize
for i in range(Norder):
    peakval = max(outlist[i, :])
    outlist[i, :] = outlist[i, :] / peakval

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)
for i in range(Norder):
    ax.plot(thetalist, outlist[i, :], label='Order'+str(orderlist[i]))
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rticks([0, 1])
ax.set_xticks(np.linspace(0, np.pi, 5))

# little bit dangerous - if filename extension is not 3 digits, there would be problem
filename = dataname[:-4] 
for i in range(Norder):
    filename = filename + '_' + str(orderlist[i])
filename = filename + '.png'
plt.legend()
plt.savefig(filename)
plt.show()'''
