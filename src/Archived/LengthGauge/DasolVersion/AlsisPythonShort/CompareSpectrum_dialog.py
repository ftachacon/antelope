#!/usr/bin/env python3
# Note that this script is only valid with display environment
# Can't be run in cluster

import sys
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

set_of_params   = sys.argv
filename = set_of_params[1]

root = tk.Tk()
root.withdraw()

Ndata = int(input('put number of files: '))
file_path = list()
for i in range(Ndata):
    file_path.append(filedialog.askopenfilename(parent=root, title='Choose a file in order you wants'))
if (Ndata != len(file_path)):
    print("Number of files is not matched - any mistake?")

print("put labels for legend")
print(len(file_path))
label_list = list()
for i in range(len(file_path)):
    label_list.append(input('label for ' +file_path[i]+': '))

width = 11
hight = width/1.62


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

for i in range(len(file_path)):
    dat = np.loadtxt(file_path[i])
    plt.plot( dat[1:, 0], np.log10(dat[1:, 5]), label=label_list[i])
plt.legend()

xp = 1.3
yp = -1.1


plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 )
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 )
plt.tick_params(labelsize = 28 )


hzero           = .98e-25

xaxmin      = np.log10(hzero)      # controling harmonic-order
xaxmax      = +5                   # controling harmonic-order

plt.ylim(xaxmin, xaxmax)
xticks0  = np.arange(1,50,4)
plt.xticks( xticks0 )

yticks0  = [-20., -15, -10.0, -5, 0, 5 ]#10., 15, 5.0];
plt.yticks( yticks0 )

plt.grid(True)


xaxmin      = 0    # controling harmonic-order axis limits, down
xaxmax      = 37   # controling harmonic-order axis limits, u


plt.xlim(xaxmin, xaxmax)
plt.ylim( np.log10(hzero),1 )

plt.savefig(filename, dpi = 300)
