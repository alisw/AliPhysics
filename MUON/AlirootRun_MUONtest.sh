#!/bin/sh
# $Id$

# first declare default values

SIMULATION=1 # will perform simulation
RECONSTRUCTION=1 # will perform reconstruction
CHECKS=1 # will perform checks
SLASHTMP=1 #will use /tmp to put the temporary raw data 
NEVENTS=100 # will simulate 100 events

#RECOPTIONS="SAVEDIGITS NOFASTDECODERS" # reconstruction options with non-high performance decoders
RECOPTIONS="SAVEDIGITS" # default reconstruction options
SIMCONFIG="$ALICE_ROOT/MUON/Config.C" # default simulation configuration file
OUTDIR=""
CURDIR=`pwd`

#RUN=0 # run number for OCDB access
SEED=1234567 # random number generator seed
SIMDIR="generated" # sub-directory where to move simulated files prior to reco
DUMPEVENT=5 # event to be dump on files

# next try to see if there are options of this script that want to change the
# defaults
 
EXIT=0

while getopts "SRX:srxn:tg:p:d:c:" option
do
  case $option in
    R ) RECONSTRUCTION=1;;
    S ) SIMULATION=1;;
    X ) 
    CHECKS=1
    DUMPEVENT=$OPTARG
    ;;
    r ) RECONSTRUCTION=0;;
    s ) SIMULATION=0;;
    x ) CHECKS=0;;
    t ) SLASHTMP=0;;
    c ) SIMCONFIG=$OPTARG;;
    d ) OUTDIR=$OPTARG;;
    n ) NEVENTS=$OPTARG;;
    g ) SEED=$OPTARG;;
    p ) RECOPTIONS=$OPTARG;; 
    *     ) echo "Unimplemented option chosen."
    EXIT=1
    ;;
  esac
done

if [ ! -n "$OUTDIR" ]; then
  OUTDIR="$CURDIR/test_out.$NEVENTS"
fi

# look if there are some leftover options
shift $(($OPTIND - 1))

if [ $# -gt 0 ] || [ "$EXIT" -eq 1 ]; then
  echo "ERROR : extra option not recognized"
  echo "Usage: `basename $0` options (-SRXsrxn:tg:p:d:c:)"
  echo "       -S (-s) perform (or not) simulation (default is do it, i.e -S)"
  echo "       -R (-r) perform (or not) reconstruction (default is do it, i.e. -R)"
  echo "       -X event (-x) perform (or not) checks and dumps (default is do it for event $DUMPEVENT, i.e. -X $DUMPEVENT)"
  echo "       -n nevents (int) number of events to simulate (default $NEVENTS)"
  echo "       -t will use OUTDIR as a tmp directory to generate raw data  "
  echo "       -g seed (uint) seed to be used in simulation (default $SEED)"
  echo "       -p recoptions (quotified string) reconstruction options to use (default \"$RECOPTIONS\")"
  echo "       -d full path to output directory (default $OUTDIR)"
  echo "       -c full path to configuration file for simulation (default $SIMCONFIG)"
  exit 4;
fi

# printout the options
echo "sim $SIMULATION rec $RECONSTRUCTION check $CHECKS"
if [ "$SIMULATION" -eq 1 ]; then
  echo "$NEVENTS events will be simulated, using the config found at $SIMCONFIG"
fi
if [ "$RECONSTRUCTION" -eq 1 ]; then
echo "Reconstruction options to be used : $RECOPTIONS"
fi
echo "Output directory will be : $OUTDIR"

if [ "$SIMULATION" -eq 1 ]; then

  rm -fr $OUTDIR
  mkdir $OUTDIR

fi

cp $ALICE_ROOT/MUON/.rootrc $ALICE_ROOT/MUON/rootlogon.C \
  $ALICE_ROOT/MUON/runReconstruction.C $ALICE_ROOT/MUON/runSimulation.C $OUTDIR

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

if [ ! -f $ALICE_ROOT/GRP/CTP/Config/Run0_999999999_v0_s1.root ]; then

  echo "Updating GRP CTP config  ..."

  aliroot -b >& $OUTDIR/updateCDBCTPConfig.out << EOF
  .L $ALICE_ROOT/MUON/UpdateCDBCTPConfig.C+
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

  aliroot -l -b -q runSimulation.C\($SEED,$NEVENTS,\""$SIMCONFIG"\"\) >& $OUTDIR/testSim.out 
  
  echo "Moving generated files to $SIMDIR"
  mkdir $OUTDIR/$SIMDIR
  mv $OUTDIR/*QA*.root $OUTDIR/*.log $OUTDIR/$SIMDIR
  mv $OUTDIR/MUON*.root $OUTDIR/Kinematics*.root $OUTDIR/galice.root $OUTDIR/TrackRefs*.root $OUTDIR/$SIMDIR

fi

###############################################################################
# 
# Performing RECONSTRUCTION
#
###############################################################################

if [ "$RECONSTRUCTION" -eq 1 ]; then

  rm -f galice.root AliESD*.root *QA*.root

  echo "Running reconstruction  ..."

  cd $OUTDIR
  
  aliroot -l -b -q runReconstruction\.C\($SEED,\""$OUTDIR/raw.root"\",\""$RECOPTIONS"\"\) >& $OUTDIR/testReco.out

fi

###############################################################################
# 
# Performing CHECKS (and dumps)
#
###############################################################################

if [ "$CHECKS" -eq 1 ]; then

  if [ -f "$OUTDIR/$SIMDIR/galice.root" ]; then

    echo "Running efficiency  ..."

    aliroot -b >& $OUTDIR/testResults.out << EOF
    .L $ALICE_ROOT/MUON/MUONefficiency.C+
    // no argument assumes Upsilon but MUONefficiency(443) works on Jpsi
    MUONefficiency("$OUTDIR/$SIMDIR/galice.root");
    .q
EOF

  if [ -f "$OUTDIR/galice.root" ]; then

      echo "Running Trigger efficiency  ..."
      aliroot -b >& $OUTDIR/testTriggerResults.out << EOF
      .L $ALICE_ROOT/MUON/MUONTriggerEfficiency.C+
      MUONTriggerEfficiency("$OUTDIR/$SIMDIR/galice.root", "$OUTDIR/galice.root", 1);
      .q
EOF

      if [ -f "$OUTDIR/AliESDs.root" ]; then

        echo "Running check ..."
        aliroot -b >& $OUTDIR/testCheck.out << EOF
        gSystem->Load("libMUONevaluation");
        .L $ALICE_ROOT/MUON/MUONCheck.C+
        MUONCheck(0, $NEVENTS-1, "$OUTDIR/$SIMDIR/galice.root", "$OUTDIR/galice.root", "$OUTDIR/AliESDs.root"); 
        .q
EOF
      fi
    fi
  fi

  echo "Running dumps for selected event ($DUMPEVENT) ..."

  if [ -f "$OUTDIR/$SIMDIR/galice.root" ]; then
    aliroot -l -b  << EOF
    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT");
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
    man->SetDefaultStorage("local://$ALICE_ROOT");
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

echo "Finished"  
echo "... see results in $OUTDIR"

cd $CURDIR
