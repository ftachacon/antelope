#!/usr/bin/env python3

# author: Dasol Kim
# date: 2022-06-27

# Important Note: This script is specialized for Honeycomb lattice BZ we used in the calculation.
#                 Case: When Γ-K distance = klength, we used x in (-3/4*klength, 3/4*klength), y in (-sqrt(3)/2*klength, sqrt(3)/2*klength)

# recall BZ structure. 
# We have reciprocal lattice vector: (0, sqrt(3)*klength), (3/2*klength, sqrt(3)/2*klength)
# This can be translated to (0, Ny-1), ((Nx-1), (Ny-1)/2) --> 2(Nx-1) symmetric in x-axis and (Ny-1) symmetric in y -axis.
# We include both end point of the boundary in momaxis.h. If the definition is changed, Nx-1 or Ny-1 should be Nx or Ny.

# To get corresponding value in arbitrary k-point, 

# description: Load data from first Brillouin zone and extend it to wider region.
#              Interpolation is used if grid is inconsistent.
#              Visulization purpose.

import enum
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.animation import FuncAnimation
import libconf, io

# import libconf and load kx, ky
def get_nx_ny(inputfile='inputParam.cfg'):
    with io.open(inputfile) as f:
        config = libconf.load(f)

    Nx = config['calc']['Npoints'][0]
    Ny = config['calc']['Npoints'][1]

    return Nx, Ny

# FIXME: This functions is hard-coded. To write general functions we also have to modify the antelope code.
def get_gridx_gridy_klength(inputfile='inputParam.cfg'):
    with io.open(inputfile) as f:
        config = libconf.load(f)

    Nx = config['calc']['Npoints'][0]
    Ny = config['calc']['Npoints'][1]

    # we set a0 = 1.0 aumstrong, and a0 is lattice paramter. (for Haldane a0 = \sqrt(3) aumstrong then.)
    a0 = 1.0 / 0.529177210903

    klength = 4*math.pi/3/a0
    gridx = np.linspace(-klength*3/4, klength*3/4, Nx)
    gridy = np.linspace(-klength*math.sqrt(3)/2, klength*math.sqrt(3)/2, Ny)
    return gridx, gridy, klength

# get index of extended BZ and return corresponding index of original BZ
# if we are in moving frame, origianl BZ is shifted by vector potential. (K(t) = k-A(t))
def index_extended_BZ(ex, ey, Nex, Ney, Nx, Ny, NAt=[0., 0.]):
    # ex, ey: index of extended BZ
    # Nex, Ney: number of grid points in extended BZ
    # Nx, Ny: number of grid points in original BZ
    # NAt: vector potential A(t) in grid index unit. (NAt = A(t)/dk)
    # return: index of original BZ (float return is allowed)

    # set zero same as original BZ
    x = ex - (Nex-1)/2 + (Nx-1)/2 - NAt[0]
    y = ey - (Ney-1)/2 + (Ny-1)/2 - NAt[1]

    # move to x-directional doubled original BZ
    x = x % (2*(Nx-1))
    y = y % (Ny-1)

    # check if x is in extended region and transform to original region
    if x > Nx-1:
        x = x - Nx
        y = (y - (Ny-1)/2) % (Ny-1)

    return x, y

# check float is almost integer
def is_integer(x, eps=1e-6):
    return abs(x - round(x)) < eps

# perform interpolation on extended BZ with data from original BZ
def arr_in_extended_BZ(arr, Nex, Ney, Nx, Ny, NAt=[0., 0.]):
    # arr: data in original BZ
    # Nex, Ney: number of grid points in extended BZ
    # Nx, Ny: number of grid points in original BZ
    # NAt: vector potential A(t) in grid index unit. (NAt = A(t)/dk)
    # return: data in extended BZ

    # set zero same as original BZ
    arr_extended = np.zeros((Nex, Ney), dtype=arr.dtype)

    # move to x-directional doubled original BZ
    for ex in range(Nex):
        for ey in range(Ney):
            x, y = index_extended_BZ(ex, ey, Nex, Ney, Nx, Ny, NAt)
            if is_integer(x) and is_integer(y):
                arr_extended[ex, ey] = arr[int(x), int(y)]
            else:
                arr_extended[ex, ey] = (1-x%1)*(1-y%1)*arr[int(x), int(y)] + (1-x%1)*(y%1)*arr[int(x), int(y)+1] + \
                     (x%1)*(1-y%1)*arr[int(x)+1, int(y)] + (x%1)*(y%1)*arr[int(x)+1, int(y)+1]

    return arr_extended

def animate_real_heat_map(original_dat, laser, shotlog, Nex, Ney, gridx, gridy, klength):
    Nx = len(gridx)
    Ny = len(gridy)
    Nshots = int(original_dat.shape[0]/Nx/Ny)
    original_dat = np.reshape(original_dat, (Nshots, Nx, Ny))

    maxdensity = np.amax(original_dat)

    gridnewx = np.linspace(-klength*3/4 * float(Nex-1)/(Nx-1), klength*3/4 * float(Nex-1)/(Nx-1), Nex)
    gridnewy = np.linspace(-klength*math.sqrt(3)/2 * float(Ney-1)/(Ny-1), klength*math.sqrt(3)/2 * float(Ney-1)/(Ny-1), Ney)

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[10,1], height_ratios=[10,1])

    ax1 = fig.add_subplot(gs[0, 0]) #heatmap
    cax1 = fig.add_subplot(gs[0, 1]) #colorbar
    ax2 = fig.add_subplot(gs[1, :]) #laser

    # def plot_init():
    #     plt.clf()
    def plot_update(i):
        # clear each axes
        ax1.clear()
        ax2.clear()
        cax1.clear()

        # plot heatmap
        new_data = arr_in_extended_BZ(original_dat[i, :, :], Nex, Ney, Nx, Ny, NAt=shotlog[i, 4:6]) # FIXME: index of vector potential is hard-coded.
        p1 = ax1.pcolormesh(gridnewx, gridnewy, np.transpose(new_data), cmap='viridis', shading='auto',vmin=0, vmax=maxdensity)
        cb1 = plt.colorbar(p1, cax=cax1)
        ax2.plot(laser[:, 0], laser[:, 1], '-b')
        ax2.scatter(shotlog[i, 1], shotlog[i, 2], c='r', marker='x')
    ani = FuncAnimation(fig, plot_update, frames=Nshots, interval=100)
    #plt.show()
    ani.save('nc.mp4', writer='ffmpeg', progress_callback=lambda i, n: print(f'Saving frame {i} / {n}'))
def main():
    Nx, Ny = get_nx_ny()
    Nex, Ney = Nx*2, int(round(Ny*1.5))

    gridx, gridy, klength = get_gridx_gridy_klength()
    gridnewx = np.linspace(-klength*3/4 * float(Nex-1)/(Nx-1), klength*3/4 * float(Nex-1)/(Nx-1), Nex)
    gridnewy = np.linspace(-klength*math.sqrt(3)/2 * float(Ney-1)/(Ny-1), klength*math.sqrt(3)/2 * float(Ney-1)/(Ny-1), Ney)

    arrPhase = [10.0*i for i in range(26)]
    foldernames = ['Phase_' + f'{arrPhase[i]:>05.1f}' for i in range(len(arrPhase))]
    savefolder = 'SetData0'
    for (i, foldername) in enumerate(foldernames):
        laser = np.loadtxt(foldername + '/outlaserdata.dat')
        shotlog = np.loadtxt(foldername + '/snapshot_log.dat')
        # nc = np.fromfile(foldername + '/first_conduction.dat')
        # pi = np.fromfile(foldername + '/first_coherence.dat', dtype=np.complex128)
        nc = np.fromfile(foldername + '/fisrt_conduction.dat') #typo...
        pi = np.fromfile(foldername + '/fisrt_coherence.dat', dtype=np.complex128)

        # just plot last shot
        Nshots = len(shotlog[:, 0])
        nc = np.reshape(nc, (Nshots, Nx, Ny))
        pi = np.reshape(pi, (Nshots, Nx, Ny))

        # plot nc
        fig = plt.figure()
        plt.pcolormesh(gridnewx, gridnewy, np.transpose(arr_in_extended_BZ(nc[-1,:,:], Nex, Ney, Nx, Ny, NAt=shotlog[-1, 4:6])), \
            cmap='viridis', shading='auto')
        plt.title('$n_{c}(k_{x},k_{y})$')
        plt.colorbar()
        plt.savefig(savefolder + '/' +  foldername + '_nc.png')
        plt.close(fig)

        # plot pi
        fig = plt.figure()
        plt.pcolormesh(gridnewx, gridnewy, np.transpose(arr_in_extended_BZ(np.abs(pi[-1,:,:]), Nex, Ney, Nx, Ny, NAt=shotlog[-1, 4:6])),\
             cmap='viridis', shading='auto')
        plt.title('$\pi_{c}(k_{x},k_{y})$')
        plt.savefig(savefolder + '/' +  foldername + '_pi.png')
        plt.close(fig)

        # plot A(t), E(t)
        # ration between maximum E(t) and A(t)
        ratio_E0_A0 = np.amax(laser[:, 1:4])/np.amax(laser[:, 4:7])
        fig, ax = plt.subplots(constrained_layout=True)
        p1 = ax.plot(laser[:,1], laser[:,2], ':b', label='E(t)')
        p2 = ax.plot(laser[:, 4] * ratio_E0_A0, laser[:, 5] * ratio_E0_A0, '--r', label='A(t)')
        ax.set_xlabel('$E_{x}$ (a.u.)')
        ax.set_ylabel('$E_{y}$ (a.u.)')
        # secondary axes are for A(t)
        secx_ax = ax.secondary_xaxis('top', functions=(lambda x: x/ratio_E0_A0, lambda x: x*ratio_E0_A0))
        secx_ax.set_xlabel('$A_{x}$ (a.u.)')
        secy_ax = ax.secondary_yaxis('right', functions=(lambda x: x/ratio_E0_A0, lambda x: x*ratio_E0_A0))
        secy_ax.set_ylabel('$A_{y}$ (a.u.)')

        ax.xaxis.label.set_color(p1[0].get_color())
        ax.yaxis.label.set_color(p1[0].get_color())
        ax.spines['left'].set_color(p1[0].get_color())
        ax.spines['bottom'].set_color(p1[0].get_color())
        ax.tick_params(axis='x', colors=p1[0].get_color())
        ax.tick_params(axis='y', colors=p1[0].get_color())

        secx_ax.xaxis.label.set_color(p2[0].get_color())
        secy_ax.yaxis.label.set_color(p2[0].get_color())
        secx_ax.spines['top'].set_color(p2[0].get_color())
        secy_ax.spines['right'].set_color(p2[0].get_color())
        secx_ax.tick_params(axis='x', colors=p2[0].get_color())
        secy_ax.tick_params(axis='y', colors=p2[0].get_color())

        plt.savefig(savefolder + '/' +  foldername + '_field.png')
        plt.close(fig)

        print("Plotting: ", i, "/", len(foldernames) )

if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='Plot extended Brillouin zone')
    # parser.add_argument('-i', '--input', type=str, default='inputParam.cfg', help='input file')

    # args = parser.parse_args()
    main()
