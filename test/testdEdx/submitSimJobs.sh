#!/bin/sh

###############################################################################
#
# Aim is to submit simulation/reconstruction jobs for multiplicity dependence check
# For local running one can use the instructions in runSim() function.
#    submitMultiplicityScan() --> main function to submit jobs for each event with a given multiplicity
#    runSim()                 --> runs the macros for simulation and reconstruction
# To see how to run these functions look below
#
###############################################################################

# global variables to be set
DATE=`date +%Y_%m_%d_%H`
nEvents=1               #fixed to 1 to make an event by event job submission
###############################################################################
submitIonTailXTalkScan()
{

  #
  # cals the submitMultiplicityScan fuction for 
  #  --> IonTail ON/OFF, XTALK ON/OFF, in reconstruction  and MC
  ###############################################################################
  ## example  --> for 5 multiplicity bins, each most central event having 15000 track multiplicity. 
  ##              For each multiplicity bin total statistic is 75000   
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh 
    submitIonTailXTalkScan 8 1000000 50 /hera/alice/marsland/software/bin/set_private_TPCdev.sh > out.log 2>&1 &
  fi
  ###############################################################################

  # inputs
  nMultBins=$1
  maxNTracks=$2
  nEventsCentral=$3
  alirootSourceFile=$4

  if [ $# -ne 4 ]; then 
    echo "4 parameters needed --> multiplicity bins, maximum number of tracks, maximum number of central events, source file for root"
    return 1
  fi

  baseDir=$(pwd)
  testDir=$baseDir/test_$nMultBins\_$maxNTracks\_$nEventsCentral\_$DATE
  DirCheckCreate $testDir
  cd $testDir

  for ((iontail=0; iontail<2 ; iontail=iontail+1))
  do
    for ((xtalk=0; xtalk<2 ; xtalk=xtalk+1))
    do

      MultiplicityScan $nMultBins $maxNTracks $nEventsCentral $alirootSourceFile $iontail $xtalk

    done  
  done

}
###############################################################################
runSim(){

  #
  # Function which runs rec.C and sim.C for given multiplicity and event number. 
  # (if submitMultiplicityScan function is called, this parameter is always fixed to 1, i.e event by event)
  # Input parameters are
  # 1) total track multiplicity 
  # 2) number of events to be processed for given total track multiplicity (if submitMultiplicityScan() is called, it is 1)
  # 3) script to source a specific aliroot
  #
  ###############################################################################
  ## example --> 2 events with total multiplicity of 1000 tracks
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh
    runSim 100 1  0 0 /hera/alice/marsland/software/bin/set_private_TPCdev.sh > out.log 2>&1 &
  fi
  ###############################################################################

  export TestdEdxNTracks=$1
  nevents=$2
  ionTail=$3
  xTalk=$4
  alirootSourceFile=$5

  if [ $# -ne 4 ]; then
    echo "5 parameters needed --> multiplicity, number of events, source file for aliroot, iontail switch (0 or 1), xTalk switch (0 or 1)"
    return 1
  fi

  echo " Running dEdx digitzer test job" 
  echo " NEvents = $nevents" 
  echo " NTracks per event  $TestdEdxNTracks"


  ## main body of the simulation part
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  echo aliroot -b -q sim.C\($nevents,$ionTail,$xTalk\)     
  aliroot -b -q sim.C\($nevents,$ionTail,$xTalk\)            2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  echo aliroot -b -q rec.C\($ionTail\,$xTalk\) 
  aliroot -b -q rec.C\($ionTail\,$xTalk\)    2>&1 | tee rec.log    
  mv syswatch.log recwatch.log

  ## OCDB entries to be dumped in human readable format
  source $ALICE_ROOT/PWGPP/CalibMacros/AliOCDBtoolkit.sh
  ocdbMakeTable AliESDs.root "ESD" OCDBrec.list
  ocdbMakeTable galice.root MC OCDBsim.list
  ocdbFileName=$(cat OCDBrec.list | grep "TPC/Calib/RecoParam" | gawk '{print $2"/"$3}' )
  dumpObject $ocdbFileName  "object" "XML" RecoParam

  return 1;

}
###############################################################################
MultiplicityScan(){
  #
  # Here we submit the jobs for the simulation//reconstruction for one setting of IonTail and XTalk configuration
  # Parameters:
  #   1. multiplicity bins to be investigated  (default 5)
  #   2. max multiplicity for whole processing (default 75000 tracks --> 5 central PbPb event )
  #   3. number of central events to be used   (default 5)
  #   4. file to source aliroot
  # (2)/(3) should be a reasonable multiplicity estimate (e.g. 15000 tracks which is 1 central PbPb event)
  # Jobs will be submitted per event     
  #
  # For each setting new directory will be created - indicating muiltiplicity
  # dir<ntracks>/dir<eventNr>  
  #
  ###############################################################################
  ## example  --> for 5 multiplicity bins, each most central event having 15000 track multiplicity. 
  ##              For each multiplicity bin total statistic is 75000   
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh 
    MultiplicityScan 8 1000000 50 /hera/alice/marsland/software/bin/set_private_TPCdev.sh 0 1 > out.log 2>&1 &
  fi
  ###############################################################################

  # inputs
  nMultBins=$1
  maxNTracks=$2
  nEventsCentral=$3
  alirootSourceFile=$4
  ionTail=$5
  xTalk=$6

  if [ $# -ne 6 ]; then 
    echo "6 parameters needed --> multiplicity bins, maximum number of tracks, maximum number of central events, source file for root, iontail switch (0 or 1), xTalk switch (0 or 1)"
    return 1
  fi

  baseDir=$(pwd)
  testDir=$baseDir/IonTail_XTalk_$ionTail\_$xTalk
  DirCheckCreate $testDir
  cd $testDir

  # create multiplicity bins 
  multPerCentralEvent=$(echo $maxNTracks/$nEventsCentral | bc)
  echo "multiplicity per most central event is $multPerCentralEvent"
  for ((i=0; i<$nMultBins ; i=i+1))
  do

    multSteps=$(echo $multPerCentralEvent/$nMultBins | bc)
    multBin=$(echo $multPerCentralEvent - $multSteps*$i | bc)
    multBinDir=$testDir/mult_$multBin
    DirCheckCreate $multBinDir
    cd $multBinDir
    echo $multBinDir

    nEventsPerMultBin=$(echo $maxNTracks/$multBin | bc)
    echo $nEventsPerMultBin
    for ((j=1; j<$(echo $nEventsPerMultBin+1 | bc) ; j=j+1))
    do

      eventDir=$multBinDir/event_$j
      DirCheckCreate $eventDir
      cd $eventDir
      cp -r $ALICE_ROOT/test/testdEdx/GRP $ALICE_ROOT/test/testdEdx/OCDB* $ALICE_ROOT/test/testdEdx/*.* .   
      #cp $ALICE_ROOT/test/testdEdx/*.* .      

      qsub -V -cwd -l h_rt=24:0:0,h_rss=6G -P alice -b y -r y -o outSim.log -e errSim.log $baseDir/submitSimJobs.sh runSim $multBin $nEvents $4 $5 $6
 
    done  
  done

  cd $baseDir 
}
###############################################################################
DirCheckCreate()
{

  #
  # check if the directory exist. If so, delete it
  #

  dirName=$1
  if [ -d "$dirName" ]; then
    echo " $dirName already exist delete it and create a new one  "
    rm -rf $dirName 
  fi
  mkdir $dirName

}
###############################################################################
main()
{
  eval "$@"
}

main $@

