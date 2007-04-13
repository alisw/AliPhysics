#!/bin/sh
rm -rf */*.root */*.dat* */*.log */fort* */hough */hlt */raw* */*~
cd ./backgr
aliroot -b -q sim.C\(2\) 2>&1 | tee sim.log
aliroot -b -q rec.C      2>&1 | tee rec.log
chmod a-w *.root
cd ../signal
aliroot -b -q sim.C\(6\) 2>&1 | tee sim.log
aliroot -b -q rec.C      2>&1 | tee rec.log


