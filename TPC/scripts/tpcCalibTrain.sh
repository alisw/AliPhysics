#
# parameters:
# 1 - basedir
# 2 - number of chunks processed
# Example: 
# source  /usr/local/grid/AliRoot/HEAD0108/TPC/scripts/tpcCalibTrain.sh `pwd`
# source /lustre/alice/marin/soft64/AliRoot/v4-17-Rev-18/TPC/scripts/tpcCalibTrain.sh `pwd`
#  work directory for test /lustre/alice/marin/rec/testRec1


#export balice=/u/miranov/.balice
#export balice=/lustre/alice/marin/soft64/setvar0417rev18.sh 
export balice=/lustre/alice/marin/soft64/setvartrunk021209.sh 
source $balice
export aliensetup=$HOME/alienSetup.sh
source $aliensetup 


#export PASS0_DIR=/usr/local/grid/AliRoot/HEAD0108
#export PASS0_DIR=/lustre/alice/marin/soft64/AliRoot/v4-17-Rev-18
export PASS0_DIR=/lustre/alice/marin/soft64/AliRoot/trunk021209
#export PASS0_DIR=$ALICE_ROOT

echo $ALICE_ROOT

#
#
# Test setup
#
export workdir=$1
export nChunks=$2
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
chgrp -R alice $workdir
chmod -R g+rwx $workdir
chmod -R o+rx $workdir
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
$ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list "alice-t3 -c 3:00 -g $bgroupNoRefit " $nChunks
nJobsNoRefitSub=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoRefit`

echo \############################   
echo Info $nJobsNoRefitSub  submitted
echo \############################   

nJobsNoDriftRun=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoDrift`

nJobsNoRefitRun=`bjobs -W | grep submitCalibJob.sh | grep -c calibNoRefit`



export totalTime=5400
export timeSleep=60
export restTime=$totalTime

while [ $restTime -gt 0 ];do
    nJobsNoDriftRun=`bjobs -g $bgroupNoDrift  -W | grep submitCalibJob.sh | grep -c calibNoDrift`
    let ratioNoDriftRunSub=100*nJobsNoDriftRun/nJobsNoDriftSub
    if [ $nJobsNoDriftSub -eq 0 ]; then
	let ratioNoDriftRunSub=0
    fi

  #  nJobsGrNoDrift=`bjobs -g $bgroupNoDrift`

    # Finding Jobs that crashed but still in queue To be debugged, does not work properly
#    comp='Break'
#    for a in `bjobs -g /recalib/calibNoDrift| gawk '{print $1}'` ; do
#	iszombie=`bpeek  $a | grep "segmentation violation"| gawk '{print $2}'`;  
#	echo Status  $iszombie , $comp
#	if [ "x$iszombie" == "x$comp" ];  then  
#	    echo The jobs id $a needs to be killed
#	    bkill  $a 
#	fi
#    done
    

    
    echo \############################   
    echo Info $nJobsNoDriftSub calibNoDrift jobs submitted  $nJobsNoDriftRun  $nJobsGrNoDrift still running $ratioNoDriftRunSub %
    echo \############################   
    
    nJobsNoRefitRun=`bjobs -g $bgroupNoRefit  -W | grep submitCalibJob.sh | grep -c calibNoRefit`
    let ratioNoRefitRunSub=100*nJobsNoRefitRun/nJobsNoRefitSub
    if [ $nJobsNoRefitSub -eq 0 ] ; then
	let ratioNoRefitRunSub=0
    fi


 #   nJobsGrNoRefit=`bjobs -g $bgroupNoRefit`

    # Finding Jobs that crashed but still in queue To be debugged does not work properly
 #   for a in `bjobs -g /recalib/calibNoRefit| gawk '{print $1}'` ; do
#	iszombie=`bpeek  $a | grep "segmentation violation"| gawk '{print $2}'`;  
#	echo Status  $iszombie , $comp
#	if [ "x$iszombie" == "x$comp" ];  then  
#	    echo The jobs id $a needs to  be killed
	#    bkill  $a 
#	fi
#    done
    


    echo \############################   
    echo Info $nJobsNoRefitSub calibNoRefit jobs  submitted  $nJobsNoRefitRun $nJobsGrNoRefit still running $ratioNoRefitRunSub %
    echo \############################   

    echo Sleeping $timeSleep , time to go $restTime    
    let restTime=restTime-timeSleep
    

   if [ $ratioNoDriftRunSub -le  10 ]; then
	if [ $ratioNoRefitRunSub -le  10 ]; then
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
$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 3:00"  "$bgroupMgNoDrift"

cd $workdir/calibNoRefit
bgroupMgNoRefit=/mergecalib/`pwd | xargs basename`
bgadd $bgroupMgNoRefit
echo  $bgroupMgNoRefit
$ALICE_ROOT/TPC/scripts/submitMerging.sh runMissing.list "alice-t3_8h -c 3:00"  "$bgroupMgNoRefit"

nJobsNoDriftMergeSub=`bjobs -W | grep -c CalibFileMerger`
nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`


export totalTime=3600
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


    if [ $ratioNoDriftMergeRunSub -le  10 ]; then
	let restTime=0
    fi


     sleep $timeSleep ;
done


echo \############################   
echo Merging  done
echo \############################   

#nJobsNoDriftMergeSub=`bjobs -W | grep -c CalibFileMerger`
#nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`


#export totalTime=3600
#export timeSleep=60
#export restTime=$totalTime

#while [ $restTime -ge 0 ];do
#    nJobsNoDriftMergeRun=`bjobs -W | grep -c CalibFileMerger`

#    echo \############################   
#    echo Info $nJobsNoDriftMergeRun jobs for merging
#    echo \############################   


#    echo Sleeping $timeSleep , time to go $restTime    
#    let restTime=restTime-timeSleep


#    if [ $nJobsNoDriftMergeRun -ge 0 ]; then
#	echo Jobs finished
#	let restTime=0;
#    fi


#    sleep $timeSleep ;
#done


echo You are done


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
echo Step 7 running tpcCalibTrain script done
echo \############################   





#cd $workdir/calibQA
#$ALICE_ROOT/TPC/scripts/makeWorkspace.sh run.list 
#$ALICE_ROOT/TPC/scripts/resubmitMissing.sh run.list alice-t3 20
#cd $workdir/
#
#
#
