#!/usr/bin/env python3
#####################################
# Warning: E field plotted in figure is x-component of laser pulse.
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import errno

shotPath = "./ShotFigures"
paramSet = np.loadtxt("setOfparameters.dat")
Nx = int(paramSet[1, 2])
Ny = int(paramSet[1, 3])
shotNumber = int(paramSet[4, 5])
Nt = int(paramSet[0, 5])
w0 = float(paramSet[0, 2])
Ncycle = int(paramSet[0, 3])

T0 = 2*math.pi/w0

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
    shotTimeStamp = np.loadtxt("densityTimestamp.dat")
    laserData = np.loadtxt('outlaserdata.dat')

    peakNumbers = list(range(-Ncycle, Ncycle+1))
    peakPositions = [T0 * i for i in peakNumbers]

    for i in range(shotNumber):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2, 2, width_ratios=[
                              10, 1], height_ratios=[10, 1])
        ax1 = fig.add_subplot(gs[0, :])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_title('density')
        p1 = ax1.imshow(density[i, :, :], origin='lower', aspect='auto', cmap=plt.get_cmap('Greens'))
        #cax1 = fig.add_subplot(gs[0, 1])
        #cb1 = plt.colorbar(p1, cax=cax1)
        p1.set_clim(0, maxdensity)
        ax2 = fig.add_subplot(gs[1, :])
        ax2.set_xticks(peakPositions)
        ax2.set_xticklabels(peakNumbers)
        ax2.set_yticklabels([])
        #markers_on = [i]
        ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
        ax2.scatter([shotTimeStamp[i, 1]], [
                    shotTimeStamp[i, 2]], color='r', marker='o')
        # ax2.plot(shotTimeStamp[:, 1], shotTimeStamp[:, 2],
        #         '+-r', markevery=markers_on)
        plt.savefig(shotPath + "/density" + "%0.8d" %
                    shotTimeStamp[i, 0] + ".png")
        plt.close(fig)
    plt.plot(laserData[:, 0], laserData[:, 1], '-b')
    plt.scatter(shotTimeStamp[:, 1],
                shotTimeStamp[:, 2], color='r', marker='o')
    plt.savefig(shotPath + "/Epulse_timestamp.png")
    plt.clf()

    polarization = np.reshape(np.fromfile(
        "polarization.dat", dtype=np.complex128), (shotNumber, Ny, Nx))
    maxpolarization = np.amax(abs(polarization))
    for i in range(shotNumber):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2, 2, width_ratios=[
                              10, 1], height_ratios=[10, 1])
        ax1 = fig.add_subplot(gs[0, :])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_title('abs square of polarization')
        p1 = ax1.imshow(
            abs(polarization[i, :, :])**2, origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
        #cax1 = fig.add_subplot(gs[0, 1])
        #cb1 = plt.colorbar(p1, cax=cax1)
        p1.set_clim(0, maxpolarization**2)
        ax2 = fig.add_subplot(gs[1, :])
        ax2.set_xticks(peakPositions)
        ax2.set_xticklabels(peakNumbers)
        ax2.set_yticklabels([])
        ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
        ax2.scatter([shotTimeStamp[i, 1]], [
                    shotTimeStamp[i, 2]], color='r', marker='o')
        plt.savefig(shotPath + "/polarization_abs2" + "%0.8d" %
                    shotTimeStamp[i, 0] + ".png")
        plt.close(fig)
    for i in range(shotNumber):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2, 2, width_ratios=[
                              10, 1], height_ratios=[10, 1])
        ax1 = fig.add_subplot(gs[0, :])
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_title('real part of polarization')
        p1 = ax1.imshow(
            np.real(polarization[i, :, :]), origin='lower', aspect='auto', cmap=cmapFull)
        #cax1 = fig.add_subplot(gs[0, 1])
        #cb1 = plt.colorbar(p1, cax=cax1)
        p1.set_clim(-maxpolarization, maxpolarization)
        ax2 = fig.add_subplot(gs[1, :])
        ax2.set_xticks(peakPositions)
        ax2.set_xticklabels(peakNumbers)
        ax2.set_yticklabels([])
        ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
        ax2.scatter([shotTimeStamp[i, 1]], [
                    shotTimeStamp[i, 2]], color='r', marker='o')
        plt.savefig(shotPath + "/polarization_real" + "%0.8d" %
                    shotTimeStamp[i, 0] + ".png")
        plt.close(fig)
    eps = 10**-18

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[
                          10, 1], height_ratios=[10, 1])
    ax1 = fig.add_subplot(gs[0, :])
    #ax1.set_title('log10(abs square of polarization)')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    p1 = ax1.imshow(
        np.log10(abs(polarization[shotNumber-1, :, :])**2 + eps), origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
    #cax1 = fig.add_subplot(gs[0, 1])
    #cb1 = plt.colorbar(p1, cax=cax1)
    #p1.set_clim(0, maxpolarization**2)
    ax2 = fig.add_subplot(gs[1, :])
    ax2.set_xticks(peakPositions)
    ax2.set_xticklabels(peakNumbers)
    ax2.set_yticklabels([])
    ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
    ax2.scatter([shotTimeStamp[shotNumber-1, 1]], [shotTimeStamp[shotNumber-1, 2]], color='r', marker='o')
    plt.savefig(shotPath + "/log_polarization_abs2" + "%0.8d" %
                shotTimeStamp[shotNumber-1, 0] + ".png")
    plt.close(fig)
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[
                          10, 1], height_ratios=[10, 1])
    ax1 = fig.add_subplot(gs[0, :])
    #ax1.set_title('log10(abs(real part of polarization))')
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    p1 = ax1.imshow(
        np.log10(abs(np.real(polarization[shotNumber-1, :, :])) + eps), origin='lower', aspect='auto', cmap=plt.get_cmap('Blues'))
    #cax1 = fig.add_subplot(gs[0, 1])
    #cb1 = plt.colorbar(p1, cax=cax1)
    #p1.set_clim(-maxpolarization, maxpolarization)
    ax2 = fig.add_subplot(gs[1, :])
    ax2.set_xticks(peakPositions)
    ax2.set_xticklabels(peakNumbers)
    ax2.set_yticklabels([])
    ax2.plot(laserData[:, 0], laserData[:, 1], '-b')
    ax2.scatter([shotTimeStamp[shotNumber-1, 1]], [shotTimeStamp[shotNumber-1, 2]], color='r', marker='o')
    plt.savefig(shotPath + "/log_polarization_real" + "%0.8d" %
                shotTimeStamp[shotNumber-1, 0] + ".png")
    plt.close(fig)
