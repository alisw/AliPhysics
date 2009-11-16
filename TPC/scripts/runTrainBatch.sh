#######################################################################
#
# Run train commands at GSI: 
# marian.ivanov@cern.ch
#
# This is just pseudo code. Bellow you can find the sequence of steps
# to be done to run and update calibration using batch farm. 
#######################################################################


####################################################################################
# Make workspace
# This is just example (see $ALICE_ROOT/TPC/scripts/ReadmeTrain.txt)
####################################################################################
cp $ALICE_ROOT/TPC/macros/CalibrateTPC.C .
cp $ALICE_ROOT/TPC/macros/ConfigOCDB.C .
#modify ConfigOCDB.C
ln -sf ~/.balice64HEAD0108 balice.sh
#use your favourite aliroot
ln -sf $ALICE_ROOT/TPC/CalibMacros/alienSetupGSI.sh alienSetup.sh
#use your alien setup 
cp $ALICE_ROOT/TPC/scripts/submitCalibJob.sh .
cp ../lists/run.list .
cp ../lists/esd.list .


####################################################################################
# 0. Create a list for each run - the superlist are located  in the lists directory
#    and make directory structure.  
#    To be in workspace dir:    
#      0.0 run.list 
#      0.1 esd.list
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list


####################################################################################
# 1. Get list of missing 
#
$ALICE_ROOT/TPC/scripts/filterMissing.sh
####################################################################################
# 2. Filter corrupted 
#
$ALICE_ROOT/TPC/scripts/filterCorrupted.sh alice-t3_8h
cat  */esd*.txt.Bad > esdBad.txt 
cat  */esd*.txt.Good > esdGood.txt 


####################################################################################
# 3. Run calibration: 
#    You have to wait until the lists are filtered
#    Only Afterwards you should process with calibration submission
#     Submitting calibration
#
$ALICE_ROOT/TPC/scripts/submitCalib.sh run.list alice-t3 50
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run85034.list alice-t3 2
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run85034.list alice-t3 5
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run90000.list alice-t3 10


###################################################################
#
# 4. Check the error and out log
# 
find `pwd`/*/err*  > errRec.log
find `pwd`/*/out*  > outRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh

###################################################################
#
# 5. Submitting merging
#
$ALICE_ROOT/TPC/scripts/submitMerging.sh run.list alice-t3_8h

###################################################################
# 
# 6. resubmit missing if neccessary
# e.g if the lists were updated
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh  run.list alice-t3 10
#
#
#$ALICE_ROOT/TPC/scripts/resubmitMissing.sh  runLaser.list alice-t3 5










#
# Merge mag field data
#
mydir=`pwd`
mkdir mergeMag
cd  mergeMag
rm mergeTrain1.txt
rm mergeTrain2.txt
touch mergeTrain1.txt
touch mergeTrain2.txt
for adir in `cat ../../lists/runMag*s.list`; do
ls $mydir/$adir/CalibObjectsTrain1.root >> mergeTrain1.txt;
ls $mydir/$adir/CalibObjectsTrain2.root >> mergeTrain2.txt;
done;
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergeTrain1.txt\"\)
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergeTrain2.txt\"\)
cd ..
#
#
# Merge mag 0
#
mydir=`pwd`
mkdir mergeMag0
cd  mergeMag0
rm mergeTrain1.txt
rm mergeTrain2.txt
touch mergeTrain1.txt
touch mergeTrain2.txt
for adir in `cat ../../lists/runMag0.list`; do
ls $mydir/$adir/CalibObjectsTrain1.root >> mergeTrain1.txt;
ls $mydir/$adir/CalibObjectsTrain2.root >> mergeTrain2.txt;
done;
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergeTrain1.txt\"\) 
aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergeTrain2.txt\"\)
cd ..
#
