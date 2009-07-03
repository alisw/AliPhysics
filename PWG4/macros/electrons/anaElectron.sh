#!/bin/bash

# ana.sh
# 
#
# Author K. Read
#
# Even though plugin version of anaElectron.C is launched using mode 4,
# this file should use mode 3 since it is executed on each GRID cpu.
#

root -b  > anaElectron.log 2>&1 <<EOF
.L anaElectron.C
anaElectron(3,"ConfigAnalysisElectron")
EOF
