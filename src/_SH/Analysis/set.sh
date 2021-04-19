#!/bin/bash

script=AnalysisReader_PyDir.py
cat dirs | while read line
do
echo "dir= $line"

"./"${script} 045 0 ${line} ${line}

done
