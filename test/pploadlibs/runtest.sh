#!/bin/bash -l
# The settings come from ~/.bash_profile

if [ "$ALICE_TARGET" = "win32gcc" ]
    then
    REXE=root_exe.exe
else
    REXE=root.exe
fi

rm -rf *.root *.dat *.log fort* hlt hough raw* *~ GRP *.ps AliHLT*

${REXE} -b -q runsim.C      2>&1 | tee sim.log
mv syswatch.log simwatch.log
${REXE} -b -q runrec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
${REXE} -b -q aod.C 2>&1 | tee aod.log
${REXE} -b -q runcheck.C 2>&1 | tee check.log



