#!/bin/bash -l
# The settings come from ~/.bash_profile

rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log GRP *.input *.data *.ps
aliroot -b -q sim.C      2>&1 | tee sim.log
mv syswatch.log simwatch.log
aliroot -b -q rec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log



