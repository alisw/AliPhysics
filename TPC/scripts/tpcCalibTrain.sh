#
# parameters:
# 1 - basedir
# 2 - number of chunks processed
# 3 - cosmic or collisions
# Example: 
# /usr/local/grid/AliRoot/HEAD0108/TPC/scripts/tpcCalibTrain.sh `pwd` 20 0 >train.log
# source /lustre/alice/marin/soft64/AliRoot/v4-17-Rev-18/TPC/scripts/tpcCalibTrain.sh `pwd`
#  work directory for test /lustre/alice/marin/rec/testRec1



export balice=/u/miranov/.balice
#export balice=/lustre/alice/marin/soft64/setvar0417rev20.sh 
#export balice=/lustre/alice/marin/soft64/setvartrunk021209.sh 
source $balice
export aliensetup=$HOME/alienSetup.sh
source $aliensetup 

echo $ALICE_ROOT

#
#
# Test setup
#
export workdir=$1
export nChunks=$2
export isCosmic=0
if [ $# -eq 3 ]; then
  isCosmic=$3
fi
echo IsCosmic  $isCosmic 



if [ ! -n length ]; then 
  echo \############################  
  echo Directory was not specified. Exiting
  echo \############################   
  return;
fi;
if [ ! -r $workdir/lists/esd.list  ] ; then
 echo \############################   
 echo File esd list does not exist. Exiting
 echo \############################   
 return;
fi; 
if [ ! -r $workdir/lists/run.list  ] ; then
 echo \############################   
 echo File run list does not exist. Exiting
 echo \############################   
 return;
fi; 

#
# Make directories
#
cd $workdir
#chgrp -R alice $workdir
#chmod -R g+rwx $workdir
#chmod -R o+rx $workdir
mkdirhier  $workdir/calibNoDrift
mkdirhier  $workdir/calibNoRefit
mkdirhier  $workdir/calibQA
#
#modify ConfigOCDB.C
#
# copy predefined Config files 
#

cp   $ALICE_ROOT/TPC/macros/CalibrateTPC.C      calibNoDrift/CalibrateTPC.C
cat  $ALICE_ROOT/TPC/macros/CalibrateTPC.C |    grep -v AddCalibCalib\(task\) > calibNoRefit/CalibrateTPC.C
cp   $ALICE_ROOT/TPC/macros/CalibrateTPC.C      calibQA/CalibrateTPC.C

if [ $isCosmic -eq 1 ] ; then
cat  $ALICE_ROOT/TPC/macros/CalibrateTPC.C|grep -v  calibTimeGain\-\>SetIsCosmic\(kFALSE\)\; > calibNoDrift/CalibrateTPC.C

cat  $ALICE_ROOT/TPC/macros/CalibrateTPC.C |grep -v AddCalibCalib\(task\) | grep -v  calibTimeGain\-\>SetIsCosmic\(kFALSE\)\; > calibNoRefit/CalibrateTPC.C

cat   $ALICE_ROOT/TPC/macros/CalibrateTPC.C |grep -v  calibTimeGain\-\>SetIsCosmic\(kFALSE\)\; >     calibQA/CalibrateTPC.C
fi

cp   $ALICE_ROOT/TPC/macros/ConfigOCDBNoDrift.C calibNoDrift/ConfigOCDB.C
cp   $ALICE_ROOT/TPC/macros/ConfigOCDBNoRefit.C calibNoRefit/ConfigOCDB.C
cp   $ALICE_ROOT/TPC/macros/ConfigOCDBQA.C      calibQA/ConfigOCDB.C
cp   lists/*.list calibNoDrift/
cp   lists/*.list calibNoRefit/
cp   lists/*.list calibQA/
ln -sf $balice          calibNoDrift/balice.sh
ln -sf $balice          calibNoRefit/balice.sh
ln -sf $balice          calibQA/balice.sh
ln -sf $aliensetup      calibNoDrift/alienSetup.sh
ln -sf $aliensetup      calibNoRefit/alienSetup.sh
ln -sf $aliensetup      calibQA/alienSetup.sh
#  make workspaces
#

echo \##################################################   
echo Info Making calibNoDrift workspace and submit jobs
echo \##################################################   


cd $workdir/calibNoDrift
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list 
bgroupNoDrift=/recalib/`pwd | xargs basename`
bgadd $bgroupNoDrift
echo  $bgroupNoDrift
bkill -g $bgroupNoDrift -r 0
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list "alice-t3 -c 3:00 -g $bgroupNoDrift " $nChunks
nJobsNoDriftSub=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoDrift`

echo \############################   
echo Info $nJobsNoDriftSub  submitted
echo \############################   

echo \##################################################   
echo Info Making calibNoRefit workspace and submit jobs
echo \##################################################   

cd $workdir/calibNoRefit
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list 
bgroupNoRefit=/recalib/`pwd | xargs basename`
bgadd $bgroupNoRefit
echo  $bgroupNoRefit
bkill -g $bgroupNoRefit -r 0
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list "alice-t3 -c 3:00 -g $bgroupNoRefit " $nChunks
nJobsNoRefitSub=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoRefit`

echo \############################   
echo Info $nJobsNoRefitSub  submitted
echo \############################   

nJobsNoDriftRun=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoDrift`

nJobsNoRefitRun=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoRefit`



export totalTime=1800
export timeSleep=60
export restTime=$totalTime

while [ $restTime -gt 0 ];do
    nJobsNoDriftRun=`bjobs -g $bgroupNoDrift  -W | grep submitCalibJob.sh | grep -c calibNoDrift`
    let ratioNoDriftRunSub=100*nJobsNoDriftRun/nJobsNoDriftSub
    if [ $nJobsNoDriftSub -eq 0 ]; then
	let ratioNoDriftRunSub=0
    fi
    
    echo \############################   
    echo Info $nJobsNoDriftSub calibNoDrift jobs submitted  $nJobsNoDriftRun  $nJobsGrNoDrift still running $ratioNoDriftRunSub %
    echo \############################   
    
    nJobsNoRefitRun=`bjobs -g $bgroupNoRefit  -W | grep submitCalibJob.sh | grep -c calibNoRefit`
    let ratioNoRefitRunSub=100*nJobsNoRefitRun/nJobsNoRefitSub
    if [ $nJobsNoRefitSub -eq 0 ] ; then
	let ratioNoRefitRunSub=0
    fi


    


    echo \############################   
    echo Info $nJobsNoRefitSub calibNoRefit jobs  submitted  $nJobsNoRefitRun $nJobsGrNoRefit still running $ratioNoRefitRunSub %
    echo \############################   

    echo Sleeping $timeSleep , time to go $restTime    
    let restTime=restTime-timeSleep
    

   if [ $ratioNoDriftRunSub -le  4 ]; then
	if [ $ratioNoRefitRunSub -le  4 ]; then
	    let restTime=0
	fi
    fi


    sleep $timeSleep ;

done

##############################################
# Submit merging when all calib jobs are done
##############################################


cd $workdir/calibNoDrift
bgroupMgNoDrift=/mergecalib/`pwd | xargs basename`
bgadd $bgroupMgNoDrift
echo $bgroupMgNoDrift
bkill -g $bgroupMgNoDrift -r 0
$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 3:00"  "$bgroupMgNoDrift"

cd $workdir/calibNoRefit
bgroupMgNoRefit=/mergecalib/`pwd | xargs basename`
bgadd $bgroupMgNoRefit
echo  $bgroupMgNoRefit
bkill -g $bgroupMgNoRefit -r 0
$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 3:00"  "$bgroupMgNoRefit"

nJobsNoDriftMergeSub=`bjobs -W | grep -c CalibFileMerger`
nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`


export totalTime=1800
export timeSleep=60
export restTime=$totalTime

while [ $restTime -gt 0 ];do
    nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`
    let ratioNoDriftMergeRunSub=100*nJobsNoDriftMergeRun/nJobsNoDriftMergeSub
    if [ $nJobsNoDriftMergeSub -eq 0 ] ; then
	let ratioNoDriftMergeRunSub=0
    fi


    echo \############################   
    echo Info $nJobsNoDriftMergeSub jobs submitted for merging, $nJobsNoDriftMergeRun running  $ratioNoDriftMergeRunSub % left
    echo \############################   


    echo Sleeping $timeSleep , time to go $restTime    
    let restTime=restTime-timeSleep


    if [ $ratioNoDriftMergeRunSub -le  4 ]; then
	let restTime=0
    fi


     sleep $timeSleep ;
done


echo \############################   
echo Merging  done
echo \############################   

#nJobsNoDriftMergeSub=`bjobs -W | grep -c CalibFileMerger`
#nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`



#################
# Step 7
#################

cd $workdir/calibNoDrift
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh run.list alice-t3_8h

cd $workdir/calibNoRefit
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh run.list alice-t3_8h

nJobsNoDriftMergeCalibSub=`bjobs -W | grep -c CalibFileMerger`
nJobsNoDriftMergeCalibRun=`bjobs -W | grep -c CalibFileMerger`

cd $workdir/calibNoDrift
find `pwd`/* | grep err > errRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh errRec.log

cd $workdir/calibNoRefit
find `pwd`/* | grep err > errRec.log
$ALICE_ROOT/TPC/scripts/filterRecLog.sh errRec.log




echo \############################   
echo MakeOCDB 
echo \############################   

cd $workdir/calibNoDrift
test -d OCDB && mv OCDB OCDB.$(date +%y%m%d_%H%M)
mkdir OCDB

runLow=$(cat ../lists/run.list | sort | head -1)
runHig=$(cat ../lists/run.list | sort | tail -1)


runL=$(echo $runLow | sed 's|^0*||')
runH=$(echo $runHig | sed 's|^0*||')

aliroot -x -q $ALICE_ROOT/TPC/CalibMacros/MakeOCDB.C\($runL,$runH,\"mergerun.list/CalibObjectsTrain1.root\"\)
cd $workdir/calibNoRefit
aliroot -x -q $ALICE_ROOT/TPC/CalibMacros/MakeOCDB.C\($runL,$runH,\"mergerun.list/CalibObjectsTrain1.root\"\)



cd $workdir/calibQA
$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list 
bgroupQA=/recalib/`pwd | xargs basename`
bgadd $bgroupQA
echo  $bgroupQA
bkill -g $bgroupQA -r 0
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list "alice-t3 -c 3:00 -g $bgroupQA " $nChunks
nJobsQASub=`bjobs -W | grep submitCalibJob.sh | grep -c calibQA`

echo \############################   
echo Info $nJobsQASub  JobsQA submitted
echo \############################   

nJobsQARun=`bjobs -W | grep submitCalibJob.sh | grep -c calibQA`



export totalTime=1800
export timeSleep=60
export restTime=$totalTime

while [ $restTime -gt 0 ];do
    nJobsQARun=`bjobs -g $bgroupQA  -W | grep submitCalibJob.sh | grep -c calibQA`
    let ratioQARunSub=100*nJobsQARun/nJobsQASub
    if [ $nJobsQASub -eq 0 ]; then
	let ratioQARunSub=0
    fi
    
    echo \############################   
    echo Info $nJobsQASub calibQA jobs submitted  $nJobsQARun  still running $ratioQARunSub %
    echo \############################   
    


    echo Sleeping $timeSleep , time to go $restTime    
    let restTime=restTime-timeSleep
    

    if [ $ratioQARunSub -le  4 ]; then
	let restTime=0
    fi


    sleep $timeSleep ;

done

##############################################
# Submit merging when the QA jobs are done
##############################################

echo \############################   
echo Going to Merge individual directories
echo \############################   


cd $workdir/calibQA
bgroupMgQA=/mergecalibQA/`pwd | xargs basename`
bgadd $bgroupMgQA
echo $bgroupMgQA
bkill -g $bgroupQA -r 0
$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 3:00"  "$bgroupMgQA"


echo \####################################   
echo Waiting for calibQA merging to finish 
echo \####################################   


nJobsQAMergeSub=`bjobs -W | grep -c CalibFileMerger`
nJobsQAMergeRun=`bjobs -W | grep -c CalibFileMerger`


export totalTime=1800
export timeSleep=60
export restTime=$totalTime

while [ $restTime -gt 0 ];do
    nJobsQAMergeRun=`bjobs -W | grep -c CalibFileMerger`
    let ratioQAMergeRunSub=100*nJobsQAMergeRun/nJobsQAMergeSub
    if [ $nJobsQAMergeSub -eq 0 ] ; then
	let ratioQAMergeRunSub=0
    fi


    echo \############################   
    echo Info $nJobsQAMergeSub jobs submitted for merging, $nJobsQAMergeRun running  $ratioQAMergeRunSub % left
    echo \############################   


    echo Sleeping $timeSleep , time to go $restTime    
    let restTime=restTime-timeSleep


    if [ $ratioQAMergeRunSub -le  4 ]; then
	let restTime=0
    fi


     sleep $timeSleep ;
done




echo \#####################################   
echo Going to do the last Merge in calibQA 
echo \#####################################   


cd $workdir/calibQA
$ALICE_ROOT/TPC/scripts/mergeCalibRun.sh run.list alice-t3_8h





echo \#####################################   
echo Last Merge in calibQA is done 
echo \#####################################   





echo \##################################   
echo MakeOCDB in the calibQA directory
echo \##################################   
# Step to be verified

cd $workdir/calibQA
test -d OCDB && mv OCDB OCDB.$(date +%y%m%d_%H%M)
mkdir OCDB

runLow=$(cat ../lists/run.list | sort | head -1)
runHig=$(cat ../lists/run.list | sort | tail -1)


runL=$(echo $runLow | sed 's|^0*||')
runH=$(echo $runHig | sed 's|^0*||')

aliroot -x -q $ALICE_ROOT/TPC/CalibMacros/MakeOCDB.C\($runL,$runH,\"mergerun.list/CalibObjectsTrain1.root\"\)




echo \#####################################   
echo Starting validation of calibQA 
echo \#####################################   
