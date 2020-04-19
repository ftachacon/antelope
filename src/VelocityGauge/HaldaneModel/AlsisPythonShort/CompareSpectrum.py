#!/usr/bin/env python3


# example:
# commandline:
# ./CompareSpectrum.py filelist.txt output.png
# ./CompareSpectrum.py filelist.txt output.png HSG  (if spectrum is HSG)
# ./CompareSpectrum.py filelist.txt output.png HHG x (x, y, lcp, rcp)

# filelist.txt:
# 
# HHG_Mt2_10.dat   "$M/t_{2}$ = 10"
# HHG_Mt2_11.dat   "$M/t_{2}$ = 11"
 
import sys
import numpy as np
import matplotlib.pyplot as plt

set_of_params   = sys.argv
argc = len(set_of_params)
listname = set_of_params[1]
filename = set_of_params[2]

spectrumType = 'HHG'
spectrumPolar = 'Total'
if (argc > 3):
  spectrumType = set_of_params[3]
if (argc > 4):
  spectrumPolar = set_of_params[4]

# for mathcing [anything, space, things in double quotes]
file_list = np.fromregex(listname, r'"(.*?)"\s+"(.+)"' ,dtype='str')
print(file_list)

width = 11
hight = width/1.62


hzero           = .98e-25

if (spectrumType == "HSG"):
  spectrum_xmin = -5
  spectrum_xmax = 15
  spectrum_ymin = np.log10(hzero)
  spectrum_ymax = 1
  xticks0  = np.arange(-4,14,4)
  yticks0  = np.arange(-20,6,5)
else:
  spectrum_xmin = 0
  spectrum_xmax = 37
  spectrum_ymin = np.log10(hzero)
  spectrum_ymax = 1
  xticks0  = np.arange(1,50,4)
  yticks0  = np.arange(-20,6,5)


fig = plt.figure(figsize=(width,hight) )
ax1 = fig.add_axes([0.2, 0.15, 0.75, 0.75])

for i in range(len(file_list[:,0])):
    dat = np.loadtxt(file_list[i,0])
    if (spectrumPolar == 'x'):
      signal = np.log10(abs(dat[1:, 1]+1j*dat[1:, 2])**2)
    elif (spectrumPolar == 'y'):
      signal = np.log10(abs(dat[1:, 3]+1j*dat[1:, 4])**2)
    elif (spectrumPolar == 'lcp'):
      signal = np.log10( abs( (dat[1:, 1]+1j*dat[1:, 2]) 
      + 1j*(dat[1:, 3]+1j*dat[1:, 4]) )**2)
    elif (spectrumPolar == 'rcp'):
      signal = np.log10( abs( (dat[1:, 1]+1j*dat[1:, 2]) 
      - 1j*(dat[1:, 3]+1j*dat[1:, 4]) )**2)
    else:
      signal = np.log10(dat[1:,5])
    plt.plot( dat[1:, 0], signal, label=file_list[i,1])
plt.legend()

plt.xlabel(r'$\rm Harmonic-Order$', fontsize=30 )
plt.ylabel(r'$\rm Log_{10}(I_{HHG})$', fontsize=30 )
plt.tick_params(labelsize = 28 )


plt.xticks( xticks0 )
plt.yticks( yticks0 )

plt.grid(True)

plt.xlim(spectrum_xmin, spectrum_xmax)
plt.ylim(spectrum_ymin, spectrum_ymax)

plt.savefig(filename, dpi = 300)
