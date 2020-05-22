#!/usr/bin/env python3

import subprocess
import io, libconf

# files to copy
configfile = 'inputParam.cfg'
jobscript = 'slurm.sh'
exefile = 'exec_hhg_mpi'
copylist = [exefile, jobscript, configfile]

# parameters to put
#arrM = [1.0+0.2*i for i in range(16)]
arrM = [-0.2*i for i in range(1,21)]

# construct folders you want
foldername = ['M_' + f'{arrM[i]:.1f}' for i in range(len(arrM))]

for i in range(len(arrM)):
    subprocess.run( ['mkdir', foldername[i]] )
    for item in copylist:
        subprocess.run( ['cp', item, foldername[i] ] , check=True)
    f = io.open(foldername[i] + '/' + configfile, 'r')
    config = libconf.load(f)
    f.close()
    config['WilsonMass']['mu'] = arrM[i]*2*config['WilsonMass']['t']
    f = io.open(foldername[i] + '/' + configfile, 'w')
    libconf.dump(config, f)

for i in range(len(arrM)):
    subprocess.run( ['sbatch', jobscript] , cwd=foldername[i])
