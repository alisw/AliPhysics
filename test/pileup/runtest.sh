#!/bin/bash -l
# The settings come from ~/.bash_profile

rm -rf *.root *.dat *.log fort* hlt hough raw* *.inp *.o GRP *.ps AliHLT*
ln -s -f ../../DPMJET/nuclear.bin
aliroot -b -q sim.C      2>&1 | tee sim.log
mv syswatch.log simwatch.log
aliroot -b -q rec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
aliroot -b -q check.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log




