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
cd $workdir
source guiEnv.sh
#
aliroot -q -b $SCRIPTDIR/ConfigOCDB.C  $ALICE_ROOT/TPC/CalibMacros/CalibEnv.C+\(\"$runList\",$startRun,$endRun\)
cp dcsTime.root $GUI_OUTDIR/time/calibTreeTime_$startRun_$endRun.root
cd $backupdir
