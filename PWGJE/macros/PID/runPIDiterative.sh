#!/bin/bash

# 1st parameter e.g. = "/lustre/nyx/alice/bhess/analysis/10d_e.pass2_merged/bhess_PID_Jets.root"
# 2nd parameter is charge mode: 0 = all, -1 = neg, +1 = pos
# 3rd parameter is TOF patching: 0 = off, 1 = on
# 4th parameter is MC: 0 = data, 1 = MC
# 5th parameter is the low entrality
# 6th parameter is the high centrality

# Modes: 0=pt, 1=z, 2=xi, 3=R, 4=jT

usePID=kTRUE  # kFALSE for prePID

for mode in 0 #1 2 3 4
do
  for jetPtRange in "-1, -1" #"10, 15" #"15, 20" "20, 30" #"10, 40" "10, 30"
  do
    aliroot runPIDiterative.C"(\"$1\", 1.8, 0.15, 50., $4, 2, 0, $usePID, $4, $mode, $2, $5, $6, $jetPtRange, 1, 1, \"\", kTRUE, kFALSE, 1, 1.0, $3)" -b -q
  done
done

# ATTENTION: For inclusive case, mode MUST bet set to 0 and $jetPtRange to -1, -1
