
from matplotlib import pyplot as plt
import numpy as np
import sys
sys.path.append('../../src/WannierGauge/LengthGauge/MovingFrame/WDSMmodel')
import postproc as pp

ws = pp.postprocessing()

current = np.loadtxt('interband_dipole_full_evol.dat').transpose()
current_local = np.loadtxt('intraband_current_full_evol.dat').transpose()

ws.plot_t_data(data_list=[current[0],current_local[0]], label_list=['$J_x$','$J_x-$'])
plt.show()

spectrum = ws.current_processing(current)
spectrum_local = ws.current_processing(current_local)

ws.plot_w_data(data_list=[spectrum,spectrum_local],label_list=['total','area'])
plt.xlim(0,20)

plt.show()
