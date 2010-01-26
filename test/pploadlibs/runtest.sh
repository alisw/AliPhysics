#!/bin/sh

if [ "$ALICE_TARGET" = "win32gcc" ]
    then
    REXE=root_exe.exe
else
    REXE=root.exe
fi

rm -rf *.root *.dat *.log fort* hlt hough raw* *~

${REXE} -b -q sim.C      2>&1 | tee sim.log
mv syswatch.log simwatch.log
${REXE} -b -q rec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
${REXE} -b -q ${ALICE_ROOT}/STEER/CreateAODfromESD.C 2>&1 | tee aod.log




