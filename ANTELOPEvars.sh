#!/bin/bash

export ANTELOPE="`dirname ${BASH_SOURCE[0]}`"
export PATH=$PATH:$ANTELOPE/utils
export PATH=$PATH:$ANTELOPE/examples
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ANTELOPE/src

