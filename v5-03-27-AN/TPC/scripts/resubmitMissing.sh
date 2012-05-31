# marian.ivanov@cern.ch
# Paramters:
# 1  - runlist
# 2  - batch queues
# 3  - number of chunks per calibration job
# 0. Find runs with esd but missing calibration
# Example:
# $ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list alice-t3 10

runlist=$1
bqueue=$2
nchunks=$3
echo runlist"      "$runlist
echo bqueue"       "$bqueue
echo nchunks"      "$nchunks

wdir=`pwd`
rm runMissing.list
touch runMissing.list
for adir in `cat $runlist`; do
  cd $wdir/$adir  
  nesd=`cat $wdir/esd$adir.txt| grep -c root`
  if [ $nesd -gt 0 ] ; then  
     ncalib=`find $wdir/$adir/ | grep -c CalibObjects`
     if [ $ncalib -lt 1 ] ; then
        echo Missing $adir
        echo $adir >> $wdir/runMissing.list
     fi;
  fi;
  cd $wdir
done;
#
# 1. Delete the content of directory
#
wdir=`pwd`
for a in `cat runMissing.list`; do
  rm -rf  $a;
done;
#
# 2. Redo directory structure for missing
#
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh runMissing.list 

#
# 3. Sumbmit calibration jobs for mssing directories
#
$ALICE_ROOT/TPC/scripts/submitCalib.sh runMissing.list "$bqueue"  $nchunks

