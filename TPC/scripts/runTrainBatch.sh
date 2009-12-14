#######################################################################
#
# Run train commands locally: 
# marian.ivanov@cern.ch
#
# This is just pseudo code. Bellow you can find the sequence of steps
# to be done to run and update calibration using batch farm. 
#######################################################################


########################################################################
# Make workspace
# This is just example (see $ALICE_ROOT/TPC/scripts/ReadmeTrain.txt)
########################################################################
cp $ALICE_ROOT/TPC/macros/CalibrateTPC.C .
cp $ALICE_ROOT/TPC/macros/ConfigOCDB.C .
#modify ConfigOCDB.C
ln -sf ~/.balice64HEAD0108 balice.sh
#use your favourite aliroot
ln -sf $HOME/alienSetup.sh alienSetup.sh
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


###########################################################################
# 1. Get list of missing 
#
##
$ALICE_ROOT/TPC/scripts/filterMissing.sh
###########################################################################
# 2. Filter corrupted 
#
$ALICE_ROOT/TPC/scripts/filterCorrupted.sh alice-t3_8h
cat  */esd*.txt.Bad > esdBad.txt 
cat  */esd*.txt.Good > esdGood.txt 


###########################################################################
# 3. Run calibration: 
#    You have to wait until the lists are filtered
#    Only Afterwards you should process with calibration submission
#     Submitting calibration
################################################################## 
################################################################## 
bgroup=/recalib/`pwd | xargs basename`
bgadd $bgroup
$ALICE_ROOT/TPC/scripts/submitCalib.sh run.list "alice-t3 -g $bgroup -c 3:00 "  20
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run85034.list alice-t3 2
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run85034.list alice-t3 5
#$ALICE_ROOT/TPC/scripts/submitCalib.sh run90000.list alice-t3 10
#$ALICE_ROOT/TPC/scripts/submitCalib.sh runMag05.list alice-t3 5
#$ALICE_ROOT/TPC/scripts/submitCalib.sh runMag02.list alice-t3 5
###################################################################
#
# 4. Check the error and out log
# 
find `pwd`/*/err*  > errRec.log
#find `pwd`/*/out*  > outRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh

###################################################################
#
# 5. Submitting merging
#
################################################################## 
################################################################## 
bgroup=/merge/`pwd | xargs basename`
bgadd $bgroup
$ALICE_ROOT/TPC/scripts/submitMerging.sh run.list "alice-t3_8h -c 0:10" $bgroup
#$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 0:10" $bgroup

###################################################################
# 
# 6. resubmit missing if neccessary
# e.g if the lists were updated
# submit in groups
# time restriction 3 hours
################################################################## 
bgroup=/recalib/`pwd | xargs basename`
bgadd $bgroup
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh  run.list "alice-t3 -c 3:00  -g $bgroup"  5
#
#
#$ALICE_ROOT/TPC/scripts/resubmitMissing.sh  runLaser.list alice-t3 5


###################################################################
# 
# 7. Merge separatelly sub run list
################################################################## 

$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh run.list 
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh runMag05.list 
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh runMag02.list
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh runMag0.list
   

ls -d mergerunMag0*.list > runMagAll.list
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh runMagAll.list
rm runMagN0.list
echo mergerunMag02.list >runMagN0.list
echo mergerunMag05.list >>runMagN0.list 
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh runMagN0.list


#
# filter debug streamers
#
rlist=runMag02.list
rm debug$rlist
for a in `cat $rlist`; do 
    ls `pwd`/$a/*/*.root  >> debug$rlist
done 
#
#
#
ls | grep Run|  sed s_Run__| sed s/_/\ /| gawk ' { print $1} '   
