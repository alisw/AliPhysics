#!/bin/sh
# $Id$

# first declare default values

SIMULATION=1 # will perform simulation
RECONSTRUCTION=1 # will perform reconstruction
RAW=1 # will reconstruct from raw data
CHECKS=1 # will perform checks
SLASHTMP=1 #will use /tmp to put the temporary raw data 
NEVENTS=50 # will simulate 100 events

#RECOPTIONS="SAVEDIGITS NOFASTDECODERS" # reconstruction options with non-high performance decoders
RECOPTIONS="SAVEDIGITS" # default reconstruction options
MC=""    # G3 Simulation with old Config.C 
#MC="g3"  # G3 Simulation (with new config macros)
#MC="g4"  # G4 Simulation (with new config macros)

OUTDIR=""
CURDIR=`pwd`

#RUN=0 # run number for OCDB access
SEED=1234567 # random number generator seed
SIMDIR="generated" # sub-directory where to move simulated files prior to reco
DUMPEVENT=5 # event to be dump on files (set to negative to skip dumps)

SIMCONFIG="$ALICE_ROOT/MUON/"$MC"Config.C"
EMBEDWITH="" # no embedding by default
REALISTIC=0 # ideal simulation by default

# next try to see if there are options of this script that want to change the
# defaults
 
EXIT=0

while getopts "SRZX:srxzn:tg:p:d:c:e:b:" option
do
  case $option in
    R ) RECONSTRUCTION=1;;
    S ) SIMULATION=1;;
    X ) 
    CHECKS=1
    DUMPEVENT=$OPTARG
    ;;
    Z ) RAW=1;;
    r ) RECONSTRUCTION=0;;
    s ) SIMULATION=0;;
    x ) CHECKS=0;;
    t ) SLASHTMP=0;;
    z ) RAW=0;;
    c ) SIMCONFIG=$OPTARG;;
    d ) OUTDIR=$OPTARG;;
    n ) NEVENTS=$OPTARG;;
    g ) SEED=$OPTARG;;
    p ) RECOPTIONS=$OPTARG;; 
    e ) EMBEDWITH=$OPTARG;;
    b ) 
    REALISTIC=$OPTARG
    RAW=0
    ;;
    *     ) echo "Unimplemented option chosen."
    EXIT=1
    ;;
  esac
done

if [ ! -n "$OUTDIR" ]; then
  if [ "$MC" = "" ]; then
    OUTDIR=$CURDIR"/test_out."$NEVENTS
  else  
    OUTDIR=$CURDIR"/"$MC"_test_out."$NEVENTS
  fi  
fi

# look if there are some leftover options
shift $(($OPTIND - 1))

if [ $# -gt 0 ] || [ "$EXIT" -eq 1 ]; then
  echo "ERROR : extra option not recognized"
  echo "Usage: `basename $0` options (-SRXsrxn:tg:p:d:c:)"
  echo "       -S (-s) perform (or not) simulation (default is do it, i.e -S)"
  echo "       -R (-r) perform (or not) reconstruction (default is do it, i.e. -R)"
  echo "       -X event (-x) perform (or not) checks and dumps (default is do it for event $DUMPEVENT, i.e. -X $DUMPEVENT)"
  echo "       -Z (-z) perform reconstruction from raw data (from digits) (default is from raw data, i.e. -Z)"
  echo "       -n nevents (int) number of events to simulate (default $NEVENTS)"
  echo "       -t will use OUTDIR as a tmp directory to generate raw data  "
  echo "       -g seed (uint) seed to be used in simulation (default $SEED)"
  echo "       -p recoptions (quotified string) reconstruction options to use (default \"$RECOPTIONS\")"
  echo "       -d full path to output directory (default $OUTDIR)"
  echo "       -c full path to configuration file for simulation (default $SIMCONFIG)"
  echo "       -e full path to a galice.root file relating to SDigits to be merged (embedding)"
  echo "       -b runnumber (int) make a realistic simulation using runnumber as anchor (default 0=ideal simulation)"
  exit 4;
fi


# printout the options
echo "sim $SIMULATION rec $RECONSTRUCTION check $CHECKS"
if [ "$SIMULATION" -eq 1 ]; then
  echo "$NEVENTS events will be simulated, using the config found at $SIMCONFIG"
fi
if [ -n "$EMBEDWITH" ]; then
  echo "Will embed simulation with $EMBEDWITH"
fi
if [ "$REALISTIC" -gt 0 ]; then
  echo "Will use anchor run $REALISTIC"
fi
if [ "$RECONSTRUCTION" -eq 1 ]; then
echo "Reconstruction options to be used : $RECOPTIONS"
if [ "$RAW" -eq 0 ]; then
echo "Will reconstruct from digits only (not from raw data)"
fi
fi
echo "Output directory will be : $OUTDIR"

if [ "$SIMULATION" -eq 1 ]; then

  rm -fr $OUTDIR
  mkdir $OUTDIR

fi

# Copy *ALL* the macros we need in the output directory, not to mess
# with our source dir in any way.
cp $ALICE_ROOT/MUON/.rootrc \
  $ALICE_ROOT/MUON/rootlogon.C \
  $ALICE_ROOT/MUON/runReconstruction.C $ALICE_ROOT/MUON/runSimulation.C \
  $ALICE_ROOT/MUON/UpdateCDBCTPConfig.C \
  $ALICE_ROOT/MUON/MUONefficiency.C \
  $ALICE_ROOT/MUON/MUONTriggerEfficiency.C \
  $ALICE_ROOT/MUON/MUONCheck.C \
  $OUTDIR

cd $OUTDIR

if [ "$SLASHTMP" -eq 0 ]; then
  mkdir ./tmp
  mkdir ./tmp/mdc1
  mkdir ./tmp/mdc2
  mkdir ./tmp/mdc1/tags

  chmod 777 ./tmp
  chmod 777 ./tmp/mdc1
  chmod 777 ./tmp/mdc2
  chmod 777 ./tmp/mdc1/tags

  export ALIMDC_RAWDB1=./tmp/mdc1
  export ALIMDC_RAWDB2=./tmp/mdc2
  export ALIMDC_TAGDB=./tmp/mdc1/tags
fi

###############################################################################
# 
# Update CTP in OCDB for MUON Trigger
#
###############################################################################

if [ ! -f $ALICE_ROOT/OCDB/GRP/CTP/Config/Run0_999999999_v0_s1.root ]; then

  echo "Updating GRP CTP config  ..."

  aliroot -b > $OUTDIR/updateCDBCTPConfig.out 2>&1 << EOF
  .L UpdateCDBCTPConfig.C++g
  UpdateCDBCTPConfig();
  .q
EOF
  
fi


###############################################################################
# 
# Performing SIMULATION
#
###############################################################################

if [ "$SIMULATION" -eq 1 ]; then

  echo "Running simulation  ..."

  aliroot -l -b -q runSimulation.C\($SEED,$NEVENTS,\""$SIMCONFIG"\"\,\""$EMBEDWITH"\"\,$REALISTIC\) > $OUTDIR/testSim.out 2>&1

  mkdir $OUTDIR/$SIMDIR

  if [ "$RAW" -eq 1 ]; then
    if [ "$REALISTIC" -eq 0 ]; then # we can not move for realistic simulations as we need e.g. kinematics to propagate the simulated vertex to the reco.
      echo "Moving generated files to $SIMDIR"
      mv $OUTDIR/*QA*.root $OUTDIR/*.log $OUTDIR/$SIMDIR
      mv $OUTDIR/MUON*.root $OUTDIR/TrackRefs*.root $OUTDIR/$SIMDIR
      mv $OUTDIR/Kinematics*.root $OUTDIR/galice.root $OUTDIR/$SIMDIR
    fi
  else  
    echo "Copying generated files to $SIMDIR"
    cp $OUTDIR/*QA*.root $OUTDIR/*.log $OUTDIR/$SIMDIR
    cp $OUTDIR/MUON*.root $OUTDIR/Kinematics*.root $OUTDIR/galice.root $OUTDIR/TrackRefs*.root $OUTDIR/$SIMDIR
  fi

  # save geometry file in a separate directory
  if [ "$MC" = "g3" ]; then
    rm -fr $ALICE_ROOT/MUON/geometry
    mkdir $ALICE_ROOT/MUON/geometry
    cp $OUTDIR/geometry.root $ALICE_ROOT/MUON/geometry
  fi 

  # copy input geometry file in a current directory
  if [ "$MC" = "g4" ]; then
    cp $ALICE_ROOT/MUON/geometry/geometry.root $OUTDIR
  fi 
  
  cp $OUTDIR/geometry.root $OUTDIR/$SIMDIR/geometry.root
  
fi

###############################################################################
# 
# Performing RECONSTRUCTION
#
###############################################################################

if [ "$RECONSTRUCTION" -eq 1 ]; then

  if [ "$RAW" -eq 1 ]; then
    if [ "$REALISTIC" -eq 0 ]; then
      rm -f galice.root
    fi
  fi  

  if [ "$REALISTIC" -ne 0 ]; then
    rm -f geometry.root
  fi
  
  rm -f AliESD*.root *QA*.root
  
  echo "Running reconstruction  ..."

  cd $OUTDIR
  
  RAWOCDB=kFALSE
  
  if [ -n "$EMBEDWITH" ]; then
    RAWOCDB=kTRUE
  fi
  
  if [ "$REALISTIC" -gt 0 ]; then
    RAWOCDB=kTRUE
  fi
  
  if [ "$RAW" -eq 1 ]; then
  
    aliroot -l -b -q runReconstruction\.C\($SEED,\""$OUTDIR/raw.root"\",\""$RECOPTIONS"\",$RAWOCDB\) > $OUTDIR/testReco.out 2>&1

  else

    aliroot -l -b -q runReconstruction\.C\($SEED,\"""\",\""$RECOPTIONS"\",$RAWOCDB\) > $OUTDIR/testReco.out  2>&1
  
  fi
  
fi

###############################################################################
# 
# Performing CHECKS (and dumps)
#
###############################################################################

if [ "$CHECKS" -eq 1 ]; then

  if [ -f "$OUTDIR/$SIMDIR/galice.root" ]; then

    echo "Running efficiency  ..."

    aliroot -b > $OUTDIR/testResults.out 2>&1 << EOF
    .L MUONefficiency.C++g
    // no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
    MUONefficiency("$OUTDIR/$SIMDIR/galice.root");
    .q
EOF

  if [ -f "$OUTDIR/galice.root" ]; then

      echo "Running Trigger efficiency  ..."
      aliroot -b > $OUTDIR/testTriggerResults.out 2>&1 << EOF
      .L MUONTriggerEfficiency.C++g
      MUONTriggerEfficiency("$OUTDIR/$SIMDIR/galice.root", "$OUTDIR/galice.root", 1);
      .q
EOF

      if [ -f "$OUTDIR/AliESDs.root" ]; then

        echo "Running check ..."
        aliroot -b > $OUTDIR/testCheck.out 2>&1 << EOF
        gSystem->Load("libMUONevaluation");
        .L MUONCheck.C++g
        MUONCheck(0, $NEVENTS-1, "$OUTDIR/$SIMDIR/galice.root", "$OUTDIR/galice.root", "$OUTDIR/AliESDs.root"); 
        .q
EOF
      fi
    fi
  fi
  
if [ "$DUMPEVENT" -ge 0 ]; then

  echo "Running dumps for selected event ($DUMPEVENT) ..."

  if [ -f "$OUTDIR/$SIMDIR/galice.root" ]; then
    aliroot -l -b  << EOF
    AliCDBManager* man = AliCDBManager::Instance();    
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliMUONMCDataInterface mcdSim("$OUTDIR/$SIMDIR/galice.root");
    mcdSim.DumpKine($DUMPEVENT);       > $OUTDIR/dump.$DUMPEVENT.kine
    mcdSim.DumpHits($DUMPEVENT);       > $OUTDIR/dump.$DUMPEVENT.hits
    mcdSim.DumpTrackRefs($DUMPEVENT);  > $OUTDIR/dump.$DUMPEVENT.trackrefs
    mcdSim.DumpDigits($DUMPEVENT,true);     > $OUTDIR/dump.$DUMPEVENT.simdigits
    mcdSim.DumpSDigits($DUMPEVENT,true);    > $OUTDIR/dump.$DUMPEVENT.sdigits
    .q
EOF
  else
    echo "$OUTDIR/$SIMDIR/galice.root is not there. Skipping sim dumps"
  fi

  if [ -f "$OUTDIR/galice.root" ]; then
    aliroot -l -b << EOF
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliMUONDataInterface dRec("$OUTDIR/galice.root");
    dRec.DumpDigits($DUMPEVENT,true); > $OUTDIR/dump.$DUMPEVENT.recdigits
    dRec.DumpRecPoints($DUMPEVENT);  > $OUTDIR/dump.$DUMPEVENT.recpoints
    dRec.DumpTrigger($DUMPEVENT); > $OUTDIR/dump.$DUMPEVENT.trigger
    .q
EOF
  else
    echo "$OUTDIR/galice.root is not there. Skipping rec dumps"
  fi
fi

fi

echo "Finished"  
echo "... see results in $OUTDIR"

cd $CURDIR
