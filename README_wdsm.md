Tested for

intel/19.1.3 impi/2019.9 mkl/2020.4

To compile wdsm version run "make" under in ./antelope; in case of success executable will be copied to ./bin/

Example input file and post processing script are in ./example/wdsm
Postprocessing is bansed on the python module (./src/WannierGauge/LengthGauge/MovingFrame/WDSMmodel/postproc.py) and helps to plot and make of fft of any 3 columns file with current evolution.

Apart from regular current/occupations files this version also produces:

- bands structure - can be plotted with postprocessing.bandstructure()
- partial currents - currents which originate from the certain area in the BZ (partial_current_w*j.dat)
- spin-currents  - 9 columns file (vxsx,vxsy,vxsz,vysx ...)
- bands resolved total current 
