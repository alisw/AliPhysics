#!/bin/sh
#
# 1 argument      - queue name
# $ALICE_ROOT/TPC/scripts/submitLaserDump.sh alice-t3_8h 

bqueue=$1
echo Queue name $bqueue
mydir=`pwd`
counterCalib=0
counterESD=0
rm runLaserESDMissing.list
for adir in `cat runLaser.list`;do
  if [ -s $mydir/$adir/esd.txt ]; then
    let counterESD=counterESD+1
  else
    echo $adir >>runLaserESDMissing.list
  fi;   
  if [ -e $mydir/$adir/CalibObjectsTrain1.root ]; then
    #ls -al $mydir/$adir/CalibObjectsTrain1.root
    cd $mydir/$adir
    bsub -q $bqueue aliroot -b -q CalibMacros/CalibLaser.C\($adir\) 
    cd $mydir 
    echo $adir
    let counterCalib=counterCalib+1
  fi;
done;
echo Laser calibration
echo Number of runs with esd= $counterESD 
echo Number of runs with calib= $counterCalib 
