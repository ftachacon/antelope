import numpy as np
import os


fname     ='BinaryTest00000_kyIndex_00000.bin';
BasicPath       = os.getcwd();

fname1         = BasicPath + fname

daty          = np.dtype("f8"); #np.dtype("i4, (3)f8")
myArray       = np.fromfile(fname1, dtype=daty);

print("\nData     =   ",myArray, "\n")
