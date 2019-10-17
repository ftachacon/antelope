#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import errno

shotPath = "./ShotFigures"
paramSet = np.loadtxt("setOfparameters.dat")
Nx = int(paramSet[1,2])
Ny = int(paramSet[1,3])
shotRate = int(paramSet[4,5])

if (shotRate > 0):
    try:
        os.mkdir(shotPath)
    #except FileExistsError:
    #    print ("Directory %s is already existed" % shotPath)
    #    pass
    except OSError as exc:
        if (exc.errno != errno.EEXIST):
            raise exc
        pass
    else:
        print ("Successfully created the directory %s " % shotPath)
    
    density = np.reshape(np.fromfile("density.dat"),(shotRate,Ny,Nx))
    shotTimeStamp = np.loadtxt("densityTimestamp.dat")
    for i in range(shotRate):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2, 2, width_ratios=[10,1], height_ratios=[10,1])
        ax1 = fig.add_subplot(gs[0,0])
        ax1.set_title('density')
        p1 = ax1.imshow(density[i, :, :], origin='lower', aspect='auto')
        cax1 = fig.add_subplot(gs[0, 1])
        cb1 = plt.colorbar(p1, cax = cax1)
        p1.set_clim(0, 0.1)
        ax2 = fig.add_subplot(gs[1,:])
        markers_on = [i]
        ax2.plot(shotTimeStamp[:, 1], shotTimeStamp[:, 2], '+-r', markevery=markers_on)
        plt.savefig(shotPath + "/density" + "%0.8d"%shotTimeStamp[i,0] + ".png")
        plt.close(fig)
    plt.plot(shotTimeStamp[:,1],shotTimeStamp[:,2], '+-r')
    plt.savefig(shotPath + "/Epulse_timestamp.png")
    plt.clf()
