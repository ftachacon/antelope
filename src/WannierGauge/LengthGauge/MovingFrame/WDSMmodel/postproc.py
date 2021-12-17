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
import io
import libconf


class postprocessing(object):
    """module for postprocessing of antelope output """

    def __init__(self, project_path='.'):
        super(postprocessing, self).__init__()
        self.project_path = project_path
        self.dim = 3
        self.width = 11
        self.hight = self.width / 1.62
        self.hzero = .98e-25
        self.spectrum_xmin = 0
        self.spectrum_xmax = 10
        self.spectrum_ymin = np.log10(self.hzero)
        self.spectrum_ymax = 1
        self.parse_input()
        self.parse_laser()
        self.color_list = ['b', 'r', 'g']
        # self.current_processing()

    def plot_t_data(self, data_list=[], label_list=[], save=False, figname=''):
        fig = plt.figure(figsize=(self.width, self.hight))
        for item, label, color in zip(data_list, label_list, self.color_list):
            plt.plot(self.t, item, lw=2., label=label, color=color)
        plt.legend(fontsize=30)
        plt.xlabel('time (a.u.) ', fontsize=30)
        plt.ylabel('$J$ (a.u.) ', fontsize=30)
        plt.tick_params(labelsize=20)
        xaxmin = self.t.min()     #
        xaxmax = self.t.max()     #
        plt.xlim(xaxmin, xaxmax)
        plt.tight_layout()
        if save:
            plt.savefig(figname)

    def fft(self, data):
        fft_current = np.zeros(data.shape, dtype=np.complex128)
        for i in range(self.dim):
            fft_current[i, :] = fft(data[i, :]) * self.dt
            fft_current[i, :] = -1j * self.w * \
                np.fft.fftshift(fft_current[i, :])
        self.fft_current = fft_current
        return self.fft_current

    def rotation(self, data):
        radiation = np.zeros(data.shape, dtype=np.complex128)
        thetaz = self.config['laser']['pulses'][0]['thetaz']
        phiz = self.config['laser']['pulses'][0]['phiz']
        radiation[0, :] = (np.cos(thetaz) * np.cos(phiz) * data[0, :]
                           + np.cos(thetaz) * np.sin(phiz) * data[1, :] - np.sin(thetaz) * data[2, :])
        radiation[1, :] = (-np.sin(phiz) * data[0, :]
                           + np.cos(phiz) * data[1, :])
        radiation[2, :] = (np.sin(thetaz) * np.cos(phiz) * data[0, :]
                           + np.sin(thetaz) * np.sin(phiz) * data[1, :]
                           + np.cos(thetaz) * data[2, :])
        self.radiation = radiation
        self.spectrum = abs(
            self.radiation[0, :])**2 + abs(self.radiation[1, :])**2 + abs(self.radiation[2, :])**2 + self.hzero
        self.spectrum_x = abs(self.radiation[0, :])**2 + self.hzero
        self.spectrum_y = abs(self.radiation[1, :])**2 + self.hzero
        self.spectrum_z = abs(self.radiation[2, :])**2 + self.hzero

    def current_processing(self, data):
        data_masked = np.zeros(data.shape, dtype=np.float64)
        for i in range(self.dim):
            data_masked[i, :] = np.blackman(self.Nt) * data[i, :]
        self.rotation(self.fft(data_masked))
        return self.spectrum, self.spectrum_x, self.spectrum_y

    def plot_w_data(self, data_list=[], label_list=[], save=False, figname=''):
        fig = plt.figure(figsize=(self.width, self.hight))
        for item, label, color in zip(data_list, label_list, self.color_list):
            # item = np.log10(item)
            plt.plot(self.w / self.w0, item, lw=2., label=label, color=color)
            # plt.fill(self.w / self.w0, item, alpha=0.1, color=color)
            plt.ylim(min(item), 1)
        plt.yscale('log')
        plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30)
        plt.ylabel(r'$\rm I_{HHG}, arb. units$', fontsize=30)
        plt.tick_params(labelsize=20)
        plt.grid(True)
        plt.xlim(self.spectrum_xmin, self.spectrum_xmax)

        plt.legend(fontsize=18)

    def parse_input(self):
        with io.open(self.project_path + '/inputParam.cfg') as f:
            self.config = libconf.load(f)

        self.Npulses = len(self.config['laser']['pulses'])
        print("\n\n==============================")
        print(f"Total {self.Npulses} pulses")
        for i in range(self.Npulses):
            print('------------------------------')
            print("Pulse ", i + 1)
            self.w0 = self.config['laser']['pulses'][i]['w0']
            print(fr"E0 = {self.config['laser']['pulses'][i]['E0']} a.u. , w0 = {self.config['laser']['pulses'][i]['w0']} a.u. ,"
                  fr"ellip = {self.config['laser']['pulses'][i]['ellip']}, ncycles = {self.config['laser']['pulses'][i]['ncycles']}, "
                  fr"cep = {self.config['laser']['pulses'][i]['cep']} rad")
            print(fr"t0 = {self.config['laser']['pulses'][i]['t0']} a.u. , "
                  fr"phix = {self.config['laser']['pulses'][i]['phix']} rad, thetaz = {self.config['laser']['pulses'][i]['thetaz']} rad, "
                  fr"phiz = {self.config['laser']['pulses'][i]['phiz']} rad, envelope = {self.config['laser']['pulses'][i]['env_name']}")
        print("==============================")

    def parse_laser(self):
        self.laser = np.transpose(np.loadtxt(
            '{}/outlaserdata.dat'.format(self.project_path)))
        self.t = self.laser[0, :]
        self.Efield = self.laser[1:3, :]
        self.Nt = len(self.t)
        self.dt = self.t[1] - self.t[0]
        self.wmax = np.pi / self.dt
        self.wmin = -self.wmax
        self.dw = 2. * self.wmax / float(self.Nt)
        self.w = np.linspace(self.wmin, self.wmax - self.dw, num=self.Nt)

    def plot_laser(self, save=False, savefig=''):
        fig = plt.figure(figsize=(self.width, self.hight))
        plt.plot(self.t, self.Efield[0, :], lw=2., color='b', label='$E_x$')
        plt.plot(self.t, self.Efield[1, :], lw=2., color='g', label='$E_y$')
        plt.legend()
        plt.xlabel('time (a.u.) ', fontsize=18)
        plt.ylabel('$E$ (a.u.) ', fontsize=18)
        plt.tick_params(labelsize=18)
        xaxmin = self.t.min()     #
        xaxmax = self.t.max()     #
        plt.xlim(xaxmin, xaxmax)
        plt.tight_layout()
        if save:
            plt.savefig(figname)
