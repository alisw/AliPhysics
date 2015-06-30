#!/bin/bash

rm -rf *.root *.dat *.log fort* hlt hough raw* GRP ITS

tar xvzf itsupcdb.tar.gz

aliroot -b -q LoadLibs.C MakeITSUSimuParam.C   2>&1 | tee simuparam.log

aliroot -b -q LoadLibs.C MakeITSRecoParam.C    2>&1 | tee recoparam.log

aliroot -b -q $1sim.C      2>&1 | tee sim.log

mv syswatch.log simwatch.log

aliroot -b -q $1rec.C      2>&1 | tee rec.log

mv syswatch.log recwatch.log



