#!/bin/bash

# Aruments
# 1   -  run list
# 2   -  start run
# 3   -  end run

runList=$1
startRun=$2
endRun=$3
echo runList=$runList
echo startRun=$startRun
echo endRun=$endRun
#
workdir=${GUI_OUTDIR}/tmp/tmp${startRun}-${endRun}
backupdir=`pwd`/
mkdirhier $workdir
cp $runList $workdir
cd $workdir
source guiEnv.sh
source $ALICE_ROOT/TPC/scripts/halloWorld.sh
#
aliroot -q -b $SCRIPTDIR/ConfigOCDB.C\($2\)  $SCRIPTDIR/CalibEnv.C+\(\"$runList\",$startRun,$endRun\)
echo End of job:
echo pwd=`pwd`
echo ls=
ls -alrt
echo cp dcsTime.root $GUI_OUTDIR/time/calibTreeTime_$startRun_$endRun.root
cp dcsTime.root $GUI_OUTDIR/time/calibTreeTime_$startRun_$endRun.root
cd $backupdir
