#!/usr/bin/env python3


# example:
# commandline:
# ./CompareSpectrum.py orderlist.txt filelist.txt output.png

# orderlist.txt:
#
# 1
# 4

# filelist.txt:
# 
# HHG_e-1.0.dat 
# HHG_e1.0.dat

# elist.txt
# -1.0
# 1.0
 
import sys
import numpy as np
import matplotlib.pyplot as plt

set_of_params   = sys.argv
argc = len(set_of_params)
ordername = set_of_params[1]
listname = set_of_params[2]
ename = set_of_params[3]
filename = set_of_params[4]

# for mathcing [anything, space, things in double quotes]
order_list = np.loadtxt(ordername)
print("order: ")
print(order_list)

e_list = np.loadtxt(ename)
print("e list: ")
print(e_list)

file_list = []
with open(listname, "r") as my_file_list:
    for i in range(len(e_list)):
        file_list.append(my_file_list.readline().replace('\n', ''))
print("file list: ")
print(file_list)


printDat = np.zeros((len(order_list), len(file_list)))

eps = 0.1
for i in range(len(file_list)):
    dat = np.loadtxt(file_list[i])
    dat = dat[1:,]
    for j in range(len(order_list)):
        temp_counter = 0.
        for k in range(len(dat[:,0])):
            if (abs(dat[k, 0] - order_list[j]) < eps):
                printDat[j, i] += dat[k, 5]
                temp_counter += 1.
            else:
                if dat[k, 0] - order_list[j] > 0:
                    break
        printDat[j, i] = printDat[j, i] / temp_counter

hzero = 1.0e-16
for i in range(len(order_list)):
    tempMax = np.amax(printDat[i, :])
    for j in range(len(file_list)):
        printDat[i, j] = np.log10(printDat[i, j] / tempMax + hzero)

for i in range(len(order_list)):
    plt.plot(np.sort(e_list), printDat[i, np.argsort(e_list)], label=r'HO = ' + str(order_list[i]))

plt.legend()

plt.show()
plt.savefig(filename, dpi = 300)
