#!/usr/bin/env python3
#####################################
# Warning: E field plotted in figure is x-component of laser pulse.
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import errno

shotPath = "./ShotFigures"
paramSet = np.loadtxt("setOfparameters.dat")
Nx = int(paramSet[1, 2])
Ny = int(paramSet[1, 3])
shotNumber = int(paramSet[4, 5])
Nt = int(paramSet[0, 5])
w0 = 6.079000e-03#float(paramSet[0, 2])
Ncycle = int(paramSet[0, 3])
a0 = paramSet[1, 4]

T0 = 2*math.pi/w0
print('T0 = ', T0)

if (shotNumber > 0):
    try:
        os.mkdir(shotPath)
    # except FileExistsError:
    #    print ("Directory %s is already existed" % shotPath)
    #    pass
    except OSError as exc:
        if (exc.errno != errno.EEXIST):
            raise exc
        pass
    else:
        print("Successfully created the directory %s " % shotPath)

    cmapFull = plt.get_cmap('bwr')
    cmapUpHalf = LinearSegmentedColormap.from_list(
        'Upper half', cmapFull(np.linspace(0.5, 1, cmapFull.N//2)))
    density = np.reshape(np.fromfile("density.dat"), (shotNumber, Ny, Nx))
    maxdensity = np.amax(density)
    mindensity = np.amin(density)
    print(maxdensity, mindensity)
    shotTimeStamp = np.loadtxt("densityTimestamp.dat")
    laserData = np.loadtxt('outlaserdata.dat')

    peakNumbers = list(range(-Ncycle, Ncycle+1))
    peakPositions = [T0 * i for i in peakNumbers]
    hzero = 0.980e-5 #* maxdensity
    hmax  =  1.0e-1;
    sig=3
    pumpAmpl = np.exp(-np.power(laserData[:,0]/T0/2.,2) )*max(laserData[:, 1])
    Kpoint = np.array([4*np.pi/3 / np.sqrt(3.)/a0, 0])
    Kdpoint = np.array([Kpoint[0]*np.cos(np.pi/3), Kpoint[0]*np.sin(np.pi/3)])
    offsetpoint = np.array([0.5, 0.5])
    rectsize = 0.05



    kgrid = np.loadtxt('mgrid.dat')
    if (kgrid[0] != 1):
        print('skip is not 1, check')
    Kx = kgrid[3:Nx+3]
    Ky = kgrid[Nx+6:Nx+6 + Ny]
    ksrefactor = 0.5
    kxmin = Kx[0]*ksrefactor
    kxmax = Kx[-1]*ksrefactor
    kymin = Ky[0]*ksrefactor
    kymax = Ky[-1]*ksrefactor
    print('KyMAx = ',kymax)
    K0=0.3853414;
    M0=0.6674/2.;
    K1 = np.array([-K0,0]);
    K2 = np.array([K0,0]);
    K3 = np.array([-K0/2.,M0]);
    K4 = np.array([K0/2.,M0]);
    K5 = np.array([-K0/2.,-M0]);
    K6 = np.array([K0/2.,-M0]);
    #print(Kx)
    #print(Ky)
    #print(shotTimeStamp)
    xshift=0.125;
    textstr=r"$\rm k-space\,\,$";
    for i in range(shotNumber):#range(2):#
        #        fig = plt.figure(constrained_layout=True)
        print ('\nSnapNo. = ',i)
        width   = 11;
        hight   = width/1.62;
        fig     = plt.figure( figsize=(width,hight) );
        #ax1     = fig.add_axes([0.10, 0.125, 0.92, 0.825])
        ax1     = plt.axes([xshift, 0.20, 0.75, 0.75])
            #gs = fig.add_gridspec(2, 2, width_ratios=[
            #                       10, 1], height_ratios=[10, 1])
        #ax1 = fig.add_subplot(gs[0, :])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        #ax1 text()
        #ax1.title(kxmin*0.9, kymax*1.2, r"$\rm k-space\,\,$", fontsize=24, color="k");
        plt.gcf().text(xshift, .96, textstr, fontsize=18,color='b')
        #ax1.set_title('log10(conduction band popluation)')
        newKx = Kx + laserData[int(shotTimeStamp[i, 0]), 3]
        newKy = Ky + laserData[int(shotTimeStamp[i, 0]), 4]
        p1 = ax1.pcolormesh(newKx, newKy, np.log10(density[i, :, :] + hzero),cmap=plt.get_cmap('seismic'))
        Apoint = (laserData[int(shotTimeStamp[i, 0]), 3], laserData[int(shotTimeStamp[i, 0]), 4])
        ax1.plot(K1[0],K1[1],'ow',markersize=10)
        ax1.plot(K4[0],K4[1],'ow',markersize=10)
        ax1.plot(K6[0],K6[1],'ow',markersize=10)
        ax1.plot(K2[0],K2[1],'o',markersize=10,color=[0,1,0.2])
        ax1.plot(K3[0],K3[1],'o',markersize=10,color=[0,1,0.2])
        ax1.plot(K5[0],K5[1],'o',markersize=10,color=[0,1,0.2])

        #krect = patches.Rectangle(Kpoint - Apoint - offsetpoint*rectsize, rectsize, rectsize, fill=False, clip_on=True)
        #ax1.add_patch(krect)
        #ax1.annotate('K', xy = Kpoint - Apoint, xytext = Kpoint - Apoint + 2*offsetpoint*rectsize, arrowprops=dict(arrowstyle='fancy'), clip_on=True)
        #kdrect = patches.Rectangle(Kdpoint - Apoint - offsetpoint*rectsize, rectsize, rectsize, fill=False, clip_on=True)
        #ax1.add_patch(kdrect)
        #ax1.annotate('K\'', xy = Kdpoint - Apoint, xytext = Kdpoint - Apoint + 2*offsetpoint*rectsize, arrowprops=dict(arrowstyle='fancy'), clip_on=True)
        ax1.set_xlim([kxmin, kxmax])
        ax1.set_ylim([kymin, kymax])
        #p1 = ax1.imshow(np.log10(density[i, :, :] + hzero), origin='lower', aspect='auto', cmap=plt.get_cmap('seismic'))
        #cax1 = fig.add_subplot(gs[0, 1])
        #cb1 = plt.colorbar(p1, cax=cax1)
        p1.set_clim(np.log10(hzero), np.log10(hmax))
        #ax2 = fig.add_subplot(gs[1, :])
        #ax1     = plt.axes([0.12, 0.25, 0.75, 0.75])
        ax2     = plt.axes([xshift, 0.01, 0.75, 0.18])
        #ax2.set_xticks(peakPositions)
        #ax2.set_xticklabels(peakNumbers)
        #ax2.set_yticklabels([])
        #markers_on = [i]
        ax2.plot(laserData[:, 0], laserData[:, 1], '-r')
        #ax2.plot(laserData[:, 0], pumpAmpl, color=[0.5,0,1.])
        #ax2.plot(laserData[:, 0], -pumpAmpl, color=[0.5,0,1.])
        ax2.fill_between( laserData[:, 0],  -pumpAmpl,  pumpAmpl, color=[0.5,0,1.], alpha=0.250 )
        
        plt.axis('off')
        ax2.scatter([shotTimeStamp[i, 1]], [
                    shotTimeStamp[i, 2]], color='b', marker='o')
        # ax2.plot(shotTimeStamp[:, 1], shotTimeStamp[:, 2],
        #         '+-r', markevery=markers_on)
        ax2.set_xlim([min(laserData[:, 0])*1.01,max(laserData[:, 0])*1.01])
        plt.savefig(shotPath + "/density" + "%0.8d" %
                    shotTimeStamp[i, 0] + ".png", dpi = 175)
        plt.close(fig)
    plt.plot(laserData[:, 0], laserData[:, 1], '-b')
    #plt.tick_params(axis='x',which='both', bottom=False, top=False, labelbotttom=False)
    plt.scatter(shotTimeStamp[:, 1],
                shotTimeStamp[:, 2], color='r', marker='o')
    plt.savefig(shotPath + "/Epulse_timestamp.png")
    plt.clf()

#    polarization = np.reshape(np.fromfile(
#        "polarization.dat", dtype=np.complex128), (shotNumber, Ny, Nx))
#    maxpolarization = np.amax(abs(polarization))
#    for i in range(shotNumber):
#        fig = plt.figure(constrained_layout=True)
#        gs = fig.add_gridspec(2, 2, width_ratios=[
#                              10, 1], height_ratios=[10, 1])
#        ax1 = fig.add_subplot(gs[0, :])
#        ax1.set_xticklabels([])
#        ax1.set_yticklabels([])
#        #ax1.set_title('log10(abs of coherence)')
#        newKx = Kx + laserData[int(shotTimeStamp[i, 0]), 3]
#        newKy = Ky + laserData[int(shotTimeStamp[i, 0]), 4]
#        p1 = ax1.pcolormesh(newKx, newKy, np.log10(abs(polarization[i, :, :] + hzero)),cmap=plt.get_cmap('PiYG'))
#        ax1.set_xlim([kxmin, kxmax])
#        ax1.set_ylim([kymin, kymax])
        #Apoint = (laserData[int(shotTimeStamp[i, 0]), 3], laserData[int(shotTimeStamp[i, 0]), 4])
        
        #krect = patches.Rectangle(Kpoint - Apoint - offsetpoint*rectsize, rectsize, rectsize, fill=False, clip_on=True)
        #ax1.add_patch(krect)
        #ax1.annotate('K', xy = Kpoint - Apoint, xytext = Kpoint - Apoint + 2*offsetpoint*rectsize, arrowprops=dict(arrowstyle='fancy'), clip_on=True)
        #kdrect = patches.Rectangle(Kdpoint - Apoint - offsetpoint*rectsize, rectsize, rectsize, fill=False, clip_on=True)
        #ax1.add_patch(kdrect)
        #ax1.annotate('K\'', xy = Kdpoint - Apoint, xytext = Kdpoint - Apoint + 2*offsetpoint*rectsize, arrowprops=dict(arrowstyle='fancy'), clip_on=True)
        #p1 = ax1.imshow(
        #    abs(polarization[i, :, :])**2, origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
        #cax1 = fig.add_subplot(gs[0, 1])
        #cb1 = plt.colorbar(p1, cax=cax1)
#        p1.set_clim(np.log10(hzero), np.log10(maxpolarization+hzero))
#        ax2 = fig.add_subplot(gs[1, :])
#        ax2.set_xticks(peakPositions)
#        ax2.set_xticklabels(peakNumbers)
#        ax2.set_yticklabels([])
#        ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
#        plt.axis('off')
#        ax2.scatter([shotTimeStamp[i, 1]], [
#                    shotTimeStamp[i, 2]], color='r', marker='o')
#        plt.savefig(shotPath + "/polarization_abs2" + "%0.8d" %
#                    shotTimeStamp[i, 0] + ".png")
#        plt.close(fig)
#    for i in range(shotNumber):
#        fig = plt.figure(constrained_layout=True)
#        gs = fig.add_gridspec(2, 2, width_ratios=[
#                              10, 1], height_ratios=[10, 1])
#        ax1 = fig.add_subplot(gs[0, :])
#        ax1.set_xticklabels([])
#        ax1.set_yticklabels([])
#        ax1.set_title('real part of polarization')
#        p1 = ax1.imshow(
#            np.real(polarization[i, :, :]), origin='lower', aspect='auto', cmap=cmapFull)
#        #cax1 = fig.add_subplot(gs[0, 1])
#        #cb1 = plt.colorbar(p1, cax=cax1)
#        p1.set_clim(-maxpolarization, maxpolarization)
#        ax2 = fig.add_subplot(gs[1, :])
#        ax2.set_xticks(peakPositions)
#        ax2.set_xticklabels(peakNumbers)
#        ax2.set_yticklabels([])
#        ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
#        ax2.scatter([shotTimeStamp[i, 1]], [
#                    shotTimeStamp[i, 2]], color='r', marker='o')
#        plt.savefig(shotPath + "/polarization_real" + "%0.8d" %
#                    shotTimeStamp[i, 0] + ".png")
#        plt.close(fig)
#    eps = 10**-18
#
#    fig = plt.figure(constrained_layout=True)
#    gs = fig.add_gridspec(2, 2, width_ratios=[
#                          10, 1], height_ratios=[10, 1])
#    ax1 = fig.add_subplot(gs[0, :])
#    #ax1.set_title('log10(abs square of polarization)')
#    ax1.set_xticklabels([])
#    ax1.set_yticklabels([])
#    p1 = ax1.imshow(
#        np.log10(abs(polarization[shotNumber-1, :, :])**2 + eps), origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
#    #cax1 = fig.add_subplot(gs[0, 1])
#    #cb1 = plt.colorbar(p1, cax=cax1)
#    #p1.set_clim(0, maxpolarization**2)
#    ax2 = fig.add_subplot(gs[1, :])
#    ax2.set_xticks(peakPositions)
#    ax2.set_xticklabels(peakNumbers)
#    ax2.set_yticklabels([])
#    ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
#    ax2.scatter([shotTimeStamp[shotNumber-1, 1]], [shotTimeStamp[shotNumber-1, 2]], color='r', marker='o')
#    plt.savefig(shotPath + "/log_polarization_abs2" + "%0.8d" %
#                shotTimeStamp[shotNumber-1, 0] + ".png")
#    plt.close(fig)
#    fig = plt.figure(constrained_layout=True)
#    gs = fig.add_gridspec(2, 2, width_ratios=[
#                          10, 1], height_ratios=[10, 1])
#    ax1 = fig.add_subplot(gs[0, :])
#    #ax1.set_title('log10(abs(real part of polarization))')
#    ax1.set_xticklabels([])
#    ax1.set_yticklabels([])
#    p1 = ax1.imshow(
#        np.log10(abs(np.real(polarization[shotNumber-1, :, :])) + eps), origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
#    #cax1 = fig.add_subplot(gs[0, 1])
#    #cb1 = plt.colorbar(p1, cax=cax1)
#    #p1.set_clim(-maxpolarization, maxpolarization)
#    ax2 = fig.add_subplot(gs[1, :])
#    ax2.set_xticks(peakPositions)
#    ax2.set_xticklabels(peakNumbers)
#    ax2.set_yticklabels([])
#    ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
#    ax2.scatter([shotTimeStamp[shotNumber-1, 1]], [shotTimeStamp[shotNumber-1, 2]], color='r', marker='o')
#    plt.savefig(shotPath + "/log_polarization_real" + "%0.8d" %
#                shotTimeStamp[shotNumber-1, 0] + ".png")
#plt.close(fig)
