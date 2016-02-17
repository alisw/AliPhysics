#!/bin/bash -l
# The settings come from ~/.bash_profile

rm -rf */*.root */*.log */*.dat */GRP */*.ps */AliHLT*
cd ./gen
aliroot -b -q rungen.C\(5\) 2>&1 | tee gen.log
chmod a-w *.root
cd ../sim
aliroot -b -q sim.C\(5\) 2>&1 | tee sim.log
aliroot -b -q rec.C      2>&1 | tee rec.log
aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log


