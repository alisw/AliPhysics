#!/bin/bash
# Aruments
# 1   -  run number
# 2   -  alien prefix

run=$1
alienPrefix=$2

echo run=$run
echo alienPrefix=$alienPrefix
#
workdir=${GUI_OUTDIR}/tmp/tmp${run}
backupdir=`pwd`/
mkdirhier $workdir
cd $workdir
source guiEnv.sh
source $ALICE_ROOT/TPC/scripts/halloWorld.sh
#
cp $SCRIPTDIR/ConfigOCDB.C .
aliroot -q -b $SCRIPTDIR/CalibSummary.C\($run\)
echo End of job:
echo pwd=`pwd`
echo ls=
ls -alrt
echo cp dcsTime.root $GUI_OUTDIR/time/calibTreeTime_$run.root
cp dcsTime.root $GUI_OUTDIR/time/calibTreeTime_$run.root


echo XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
echo Copy ALIEN prefix $alienPrefix 
echo XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
echo alien_rm $alienPrefix/calibTreeTime_$run.root
echo alien_cp `pwd`/dcsTime.root $alienPrefix/calibTreeTime_$run.root
alien_rm $alienPrefix/calibTreeTime_$run.root
alien_cp dcsTime.root $alienPrefix/calibTreeTime_$run.root

cd $backupdir
