#!/bin/bash -l
# The settings come from ~/.bash_profile

rm -rf */*.root */*.dat* */*.log */fort* */hough */hlt */raw* */*~ */GRP */*.ps */AliHLT*
cd ./backgr
aliroot -b -q sim.C\(2\) 2>&1 | tee sim.log
aliroot -b -q rec.C      2>&1 | tee rec.log
aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log
chmod a-w *.root
cd ../signal
aliroot -b -q sim.C\(6\) 2>&1 | tee sim.log
aliroot -b -q rec.C      2>&1 | tee rec.log
aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
aliroot -b -q aod.C 2>&1 | tee aod.log


