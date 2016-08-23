#!/bin/bash -l
# The settings come from ~/.bash_profile

# Clean-up
rm -rf */*.root */*.log */*.dat */GRP */*.ps */AliHLT* */*.d */*.pcm */*.so

# Generation of kinematics tree
cd ./gen
root.exe -b -q rungen.C\(50\) 2>&1 | tee gen.log
chmod a-w *.root

# Geant3 simulation
cd ../sim
root.exe -b -q runsim.C\(50\) 2>&1 | tee sim.log
mv syswatch.log simwatch.log
root.exe -b -q runrec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
root.exe -b -q runcheck.C 2>&1 | tee check.log
root.exe -b -q aod.C 2>&1 | tee aod.log

# Geant4 simulation
cd ../simG4
root.exe -b -q runsim.C\(50\) 2>&1 | tee sim.log
mv syswatch.log simwatch.log
root.exe -b -q runrec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
root.exe -b -q runcheck.C 2>&1 | tee check.log
root.exe -b -q aod.C 2>&1 | tee aod.log


