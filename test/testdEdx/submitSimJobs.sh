#!/bin/bash

###############################################################################
#
# Aim is to submit simulation/reconstruction jobs for multiplicity dependence check
# Functions
#    IonTailXTalkScan()       --> runs over IonTail ON/OFF, XTALK ON/OFF cases in rec.C and sim.C 
#    submitMultiplicityScan() --> main function to submit jobs for each event with a given multiplicity
#    runSim()                 --> runs the macros for simulation and reconstruction
#    RunCallgrind()           --> runs callgrind for a given setting of IonTail ON/OFF, XTALK ON/OFF
#    RunValgrind()            --> runs valgrind for a given setting of IonTail ON/OFF, XTALK ON/OFF
# To see how to run each functions individually look at the instructions written in the function itself
#
###############################################################################

###############################################################################
#### global variables to be set
onlyValgrind=0                                                                      # to run only callgrind and valgrind
nEventsPerJob=1                                                                     # fixed to 1 to make an event by event job submission
valgrindNEvents=1                                                                   # fixed to 1 for fast output valgrind test
valgrindNTracks=100                                                                 # fixed to 50 for fast output valgrind test
#alirootSource="source /hera/alice/marsland/software/bin/set_private_TPCdev.sh"      # either an alias or a "source" statement
alirootSource=""      # either an alias or a "source" statement
callgrindCommand="/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes"
valgrindCommand="/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v"
qsubCommand="qsub -V -cwd -l h_rt=24:0:0,h_rss=6G -P alice -b y -r y -o outSim.log -e errSim.log"
workDir=`pwd`

## meaning of the lookup array settings [ON:1, OFF:0]
# 0 --> XtalkMC ON, others OFF
# 1 --> XtalkMC and XtalkRec ON, others OFF
# 2 --> IonTailMC ON, others OFF
# 3 --> IonTailMC and IonTailRec ON, others OFF
# 4 --> ALL OFF
# 5 --> All ON
arrayXTalkSwitchMC=(1 1 0 0 0 1)
arrayXTalkSwitchRec=(0 1 0 0 0 1)
arrayIonTailSwitchMC=(0 0 1 1 0 1)
arrayIonTailSwitchRec=(0 0 0 1 0 1)
DATE=`date +%Y%m%d_%H`                                                              # for naming of a specific test
###############################################################################

###############################################################################
####################### Rest of the Script is generic #########################
###############################################################################
IonTailXTalkScan()
{

  #
  # Scans over the cases of IonTail:ON/OFF, XTALK:ON/OFF, and for each calls submitMultiplicityScan()
  ###############################################################################
  # Parameters:
  # 1) number os multiplicity bins
  # 2) total track multiplicity which is fixed for each multbin
  # 3) number of central events to be analysed
  #   (e.g. (2)=75000, (3)=5 for a central event of 15000) 
  ###############################################################################
  # how to run:
  if [ 1 -eq 0 ]; then
    source submitSimJobs.sh 
    IonTailXTalkScan 2 1000 2 > outIonTailXTalkScan.log 2>&1 &
  fi
  ###############################################################################

  # inputs
  nMultBins=$1
  maxNTracks=$2
  nEventsCentral=$3

  if [ $# -ne 3 ]; then
    echo "IonTailXTalkScan: 3 parameters needed --> multiplicity bins, maximum number of tracks, maximum number of central events"
    return 1
  fi

  baseDir=$(pwd)
  testDir=$baseDir/test_$nMultBins\_$maxNTracks\_$nEventsCentral\_$DATE
  DirCheckCreate $testDir
  cd $testDir

  arrayLength=${#arrayXTalkSwitchMC[@]} 
  for ((i=0; i<$arrayLength ; i=i+1))
  do
    $baseDir/submitSimJobs.sh MultiplicityScan $nMultBins $maxNTracks $nEventsCentral ${arrayXTalkSwitchMC[$i]} ${arrayXTalkSwitchRec[$i]} ${arrayIonTailSwitchMC[$i]} ${arrayIonTailSwitchRec[$i]}
    echo " settings are = " $nMultBins $maxNTracks $nEventsCentral ${arrayXTalkSwitchMC[$i]} ${arrayXTalkSwitchRec[$i]} ${arrayIonTailSwitchMC[$i]} ${arrayIonTailSwitchRec[$i]}

  done

}
###############################################################################
MultiplicityScan(){
  #
  # Here we submit the jobs for the simulation//reconstruction for one setting of IonTail and XTalk configuration
  # Parameters:
  # 1) multiplicity bins to be investigated  (default 5)
  # 2) max multiplicity for whole processing (default 75000 tracks --> 5 central PbPb event )
  # 3) number of central events to be used   (default 5)
  # 4) Xtalk MC switch
  # 5) Xtalk Rec switch
  # 6) IonTail MC switch
  # 7) IonTail Rec switch
  # (2)/(3) should be a reasonable multiplicity estimate (e.g. 15000 tracks which is 1 central PbPb event)
  # Jobs will be submitted per event
  #
  # For each setting new directory will be created - indicating muiltiplicity
  # dir<ntracks>/dir<eventNr>
  #
  ###############################################################################
  ## how to run  --> for 5 multiplicity bins, each most central event having 15000 track multiplicity.
  ##                 For each multiplicity bin total statistic is 75000
  ##                 xtalkMC is ON, xtalkRec is OFF, iontailMC is OFF, iontailRec OFF
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh
    MultiplicityScan 2 1000000 2 1 0 0 0 > outMultiplicityScan.log 2>&1 &
  fi
  ###############################################################################

  # inputs
  nMultBins=$1
  maxNTracks=$2
  nEventsCentral=$3
  xTalkMCswitch=$4
  xTalkRecswitch=$5
  ionTailMCswitch=$6
  ionTailRecswitch=$7


  if [ $# -ne 7 ]; then
    echo "MultiplicityScan: 7 parameters needed --> multiplicity bins, maximum number of tracks, maximum number of central events, xTalkMCswitch, xTalkRecswitch, ionTailMCswitch, ionTailRecswitch"
    return 1
  fi

  baseDir=$(pwd)
  testDir=$baseDir/XTalk_mc$xTalkMCswitch\-rec$xTalkRecswitch\_IONTAIL_mc$ionTailMCswitch\-rec$ionTailRecswitch
  DirCheckCreate $testDir
  cd $testDir

  echo " ====================================================== "
  echo " first submit valgrind and callgrind jobs for debugging "
  valDir=$testDir/valgrind; DirCheckCreate $valDir; cd $valDir; cp -r $workDir/*.* .
  eval $qsubCommand $workDir/submitSimJobs.sh RunValgrind $valgrindNTracks $valgrindNEvents $4 $5 $6 $7 
  callDir=$testDir/callgrind; DirCheckCreate $callDir; cd $callDir; cp -r $workDir/*.* .
  eval $qsubCommand $workDir/submitSimJobs.sh RunCallgrind  $valgrindNTracks $valgrindNEvents $4 $5 $6 $7 
  echo " ====================================================== "

  if [ $onlyValgrind == 0 ]; then
    # create multiplicity bins
    cd $testDir
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
        cp -r $workDir/*.* .

        eval $qsubCommand $workDir/submitSimJobs.sh runSim $multBin $nEventsPerJob $4 $5 $6 $7

      done
    done
  fi

  cd $baseDir
}
###############################################################################
runSim(){

  #
  # Function which runs rec.C and sim.C for given multiplicity and event number.
  # (if submitMultiplicityScan function is called, this parameter is always fixed to 1, i.e event by event)
  # Input parameters are
  # 1) total track multiplicity 
  # 2) number of events to be processed for given total track multiplicity (if submitMultiplicityScan() is called, it is 1)
  # 3) Xtalk MC switch
  # 4) Xtalk Rec switch
  # 5) IonTail MC switch
  # 6) IonTail Rec switch
  ###############################################################################
  ## how to run --> 2 events with total multiplicity of 1000 tracks
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh
    runSim 100 1 1 0 0 0  > outRunSim.log 2>&1 &
  fi
  ###############################################################################

  export TestdEdxNTracks=$1
  nEventsPerJob=$2
  xTalkMCswitch=$3
  xTalkRecswitch=$4
  ionTailMCswitch=$5
  ionTailRecswitch=$6

  if [ $# -ne 6 ]; then
    echo "runSim: 6 parameters needed --> multiplicity, number of events, xTalkMCswitch, xTalkRecswitch, ionTailMCswitch, ionTailRecswitch"
    return 1
  fi

  echo " Running dEdx digitzer test job" 
  echo " NEvents = $nEventsPerJob"
  echo " NTracks per event  $TestdEdxNTracks"

  # source aliroot environment
  eval $alirootSource
  echo " ==================== ALIROOT environment used ======================== "
  cd $ALICE_ROOT
  which aliroot
  echo $(pwd)
  git describe --all 
  git describe --dirty
  cd -
  echo " ====================================================================== "


  ## main body of the simulation part
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  printf   "\n ======================================================================\n\n"
  echo Running: aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)        
  aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)                                      #  make a specific OCDB for simulation
  
  echo Running: aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)
  aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)           2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  printf   "\n ======================================================================\n\n"
  
  echo Running: aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)
  aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)    2>&1 | tee rec.log    
  mv syswatch.log recwatch.log
  
  ## OCDB entries to be dumped in human readable format
  source $ALICE_ROOT/PWGPP/CalibMacros/AliOCDBtoolkit.sh
  printf   "\n ======================================================================\n\n"
  echo Running: ocdbMakeTable AliESDs.root "ESD" OCDBrec.list
  ocdbMakeTable AliESDs.root "ESD" OCDBrec.list
  printf   "\n ======================================================================\n\n"
  echo Running: ocdbMakeTable galice.root MC OCDBsim.list
  ocdbMakeTable galice.root MC OCDBsim.list
  ocdbFileName=$(cat OCDBrec.list | grep "TPC/Calib/RecoParam" | gawk '{print $2"/"$3}' )
  printf   "\n ======================================================================\n\n"
  echo Running: dumpObject $ocdbFileName  "object" "XML" RecoParam
  dumpObject $ocdbFileName  "object" "XML" RecoParam

  return 1;

}
###############################################################################
RunCallgrind()
{
  # Run Callgrind
  # Input parameters are
  # 1) total track multiplicity 
  # 2) number of events to be processed for given total track multiplicity 
  # 3) Xtalk MC switch; 0 or 1
  # 4) Xtalk Rec switch; 0 or 1
  # 5) IonTail MC switch; 0 or 1
  # 6) IonTail Rec switch; 0 or 1
  ###############################################################################
  ## how to run 
  if [ 1 -eq 0 ]; then
    source submitSimJobs.sh
    RunCallgrind 100 1 1 0 0 0 > outRunCallgrind.log 2>&1 &
  fi
  ###############################################################################

  export TestdEdxNTracks=$1
  nEventsPerJob=$2
  xTalkMCswitch=$3
  xTalkRecswitch=$4
  ionTailMCswitch=$5
  ionTailRecswitch=$6

  if [ $# -ne 6 ]; then
    echo "runSim: 6 parameters needed --> multiplicity, number of events, xTalkMCswitch, xTalkRecswitch, ionTailMCswitch, ionTailRecswitch"
    return 1
  fi

  eval $alirootSource

  ## main body of the simulation part
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  printf   "\n ======================================================================\n\n"
  echo Running: aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)        
  aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)                                      #  make a specific OCDB for simulation
  
  # run callgrind for sim.C
  echo Running: callgrind aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)
  $callgrindCommand aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)            2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  printf   "\n ======================================================================\n\n"

  # new directory for valgrind of the rec
  callDirRec=$(pwd)/callgrind_Rec
  mkdir $callDirRec
  cp -rf $(pwd)/GRP $(pwd)/OCDB*  $(pwd)/*.*  $callDirRec/
  cd $callDirRec

  # run callgrind for rec.C
  echo Running: callgrind aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)
  $callgrindCommand aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)                         2>&1 | tee rec.log    
  mv syswatch.log recwatch.log

  return 1;

}
###############################################################################
###################
RunValgrind()
{
#
# Run Valgrind 
# Input parameters are
# Input parameters are
# 1) total track multiplicity 
# 2) number of events to be processed for given total track multiplicity (if submitMultiplicityScan() is called, it is 1)
# 3) Xtalk MC switch
# 4) Xtalk Rec switch
# 5) IonTail MC switch
# 6) IonTail Rec switch
###############################################################################
## how to run 
  if [ 1 -eq 0 ]; then
    cd $ALICE_ROOT/test/testdEdx
    source submitSimJobs.sh
    RunValgrind 100 1 1 0 0 0 > outRunValgrind.log 2>&1 &
  fi
###############################################################################

  export TestdEdxNTracks=$1
  nEventsPerJob=$2
  xTalkMCswitch=$3
  xTalkRecswitch=$4
  ionTailMCswitch=$5
  ionTailRecswitch=$6

  if [ $# -ne 6 ]; then
    echo "runSim: 6 parameters needed --> multiplicity, number of events, xTalkMCswitch, xTalkRecswitch, ionTailMCswitch, ionTailRecswitch"
    return 1
  fi

  eval $alirootSource

  ## main body of the simulation part
  rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
  printf   "\n ======================================================================\n\n"
  echo Running: aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)        
  aliroot -b -q sim.C\(-1,$ionTailMCswitch,$xTalkMCswitch\)                                      #  make a specific OCDB for simulation

  # run valgrind for sim.C  
  echo Running: valgrind aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)
  $valgrindCommand aliroot -b -q sim.C\($nEventsPerJob,$ionTailMCswitch,$xTalkMCswitch\)    2>&1 | tee sim.log
  mv syswatch.log simwatch.log
  printf   "\n ======================================================================\n\n"
 
  # new directory for valgrind of the rec
  valDirRec=$(pwd)/valgrind_Rec
  mkdir $valDirRec
  cp -r $(pwd)/GRP $(pwd)/OCDB*  $(pwd)/*.*  $valDirRec/
  cd $valDirRec

  # run valgrind for rec.C
  echo Running: valgrind aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)
  $valgrindCommand aliroot -b -q rec.C\($ionTailRecswitch\,$xTalkRecswitch\)                 2>&1 | tee rec.log    
  mv syswatch.log recwatch.log
  
  return 1;
  
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

