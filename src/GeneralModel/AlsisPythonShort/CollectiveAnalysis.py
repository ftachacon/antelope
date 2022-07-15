#!/usr/bin/env python3

# %%

import subprocess
from pathlib import Path

# script for analysis
#analysisScript = './AnalysisReader_PyDir.py'
analysisScript = './AnalysisInterIntra.py'
# for checking
checkFile =  'interband_dipole_full_evol.dat'
# input file
configfile = 'inputParam.cfg'

# get the list of sub-directories
currentPath = Path.cwd()
dirList = [x for x in currentPath.iterdir() if x.is_dir() and (x/checkFile).exists()]

for dr1 in currentPath.iterdir():
    if dr1.is_dir():
        print(dr1.name)
        result = subprocess.run([analysisScript, dr1.name, dr1.name, '--noshow', '--f', '-mn' ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print(result.stdout)
        print(result.stderr)
