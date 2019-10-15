#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import errno

shotPath = "./ShotFigures"
paramSet = np.loadtxt("setOfparameters.dat")
Nx = int(paramSet[1,2])
Ny = int(paramSet[1,3])
shotRate = int(paramSet[4,5]) + 1

if (shotRate > 0):
    try:
        os.mkdir(shotPath)
    #except FileExistsError:
    #    print ("Directory %s is already existed" % shotPath)
    #    pass
    except OSError as exc:
        if (exc.errno != errno.EEXIST):
            raise exc;
        pass
    else:
        print ("Successfully created the directory %s " % shotPath)
    
    density = np.reshape(np.fromfile("density.dat"),(shotRate,Ny,Nx))
    shotTimeStamp = np.loadtxt("densityTimestamp.dat")
    for i in range(shotRate):
        plt.imshow(density[i,:,:])
        plt.colorbar()
        plt.savefig(shotPath + "/density" + "%0.8d"%shotTimeStamp[i,0] + ".png")
        plt.clf()
    plt.plot(shotTimeStamp[:,1],shotTimeStamp[:,2], '+-r')
    plt.savefig(shotPath + "/Epulse_timestamp.png")
    plt.clf()
