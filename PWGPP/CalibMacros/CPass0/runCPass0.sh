#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration 
#
# Files assumed to be in working directory:
# recCPass0.C          - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments (run locally):
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 

# example:
# runCPass0.sh raw.root  50  104892

#ALIEN setting
# $1 = raw input filename
runNum=`echo $1 | cut -d "/" -f 6 | sed 's/^0*//'`

# Exporting variable to define that we are in CPass0 to be used in reconstruction
export CPass='0'

# Set memory limits to a value lower than the hard site limits to at least get the logs of failing jobs
ulimit -S -v 3500000

if [ $# -lt 3 ]; then
    # alien Setup
    nEvents=99999999
    ocdbPath="raw://"
    
    # use the value passed by LPM, or by default use the kCalibBarrel alias
    triggerAlias=${ALIEN_JDL_TRIGGERALIAS-?Trigger=kCalibBarrel}
    #triggerAlias="?Trigger=kPhysicsAll"
fi

if [ $# -ge 4 ]; then
  # local setup
  nEvents=$2
  runNum=$3
  ocdbPath=$4
  triggerAlias="?Trigger=kCalibBarrel"
fi

if [ $# -eq 5 ]; then
    # local setup in case we provide the trigger mask
    # the trigger mask is first stripped of quotation characters
    triggerAlias=${5//\"/}
fi

# ===| TPC default values |===================================================
# can be overwritten by JDL below
# JDL will overwrite the config file
#
# ---| gain calibration |-----------------------------------------------------
# default in CPass0 is full calibration,
#    for number convention see AliTPCPreprocessorOffline::EGainCalibType
#
export TPC_CPass0_GainCalibType=1

# ===| TPC JDL overwrites |===================================================
#
export TPC_CPass0_GainCalibType=${ALIEN_JDL_TPC_CPASS0_GAINCALIBTYPE-$TPC_CPass0_GainCalibType}

echo "TPC_CPass0_GainCalibType=${TPC_CPass0_GainCalibType}" | tee -a calib.log

export TPC_GainCalib_minSignalN=${ALIEN_JDL_TPC_GAINCALIB_MINSIGNALN-$TPC_GainCalib_minSignalN}

echo "TPC_GainCalib_minSignalN=${TPC_GainCalib_minSignalN}" | tee -a calib.log

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi

echo "File to be  processed $1"
echo "Number of events to be processed $nEvents"

echo "* PATH: $PATH"
echo "* LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo

if [ "$2" == "OCDB" ]; then
    echo "Generating OCDB.root only"
    export OCDB_SNAPSHOT_CREATE="kTRUE"
    export OCDB_SNAPSHOT_FILENAME="OCDB.root"
    touch OCDB.generating.job
fi

CHUNKNAME="$1"

if [ "${CHUNKNAME:0:1}" = "/" ]; then
    FILENAME=${CHUNKNAME##*/}

    if [ -f "$FILENAME" ]; then
        # locally downloaded chunk
        CHUNKNAME="`pwd`/$FILENAME"
    else
        # one chunk from alien (nodownload option to the collection)
        CHUNKNAME="alien://$CHUNKNAME"
    fi
fi

if [ -f "wn.xml" ]; then
    CHUNKNAME="collection://wn.xml"
fi

echo "* Running AliRoot to reconstruct $*"
echo "* Chunk name: $CHUNKNAME"
echo "* Run number: $runNum"
echo "* nEvents: $nEvents"
echo "* ocdbPath: $ocdbPath"
echo "* triggerAlias: $triggerAlias"
echo "* additionalRecOptions: $additionalRecOptions"
echo ""

echo executing aliroot -l -b -q -x "recCPass0.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerAlias\", \"$additionalRecOptions\")"
time aliroot -l -b -q -x "recCPass0.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerAlias\", \"$additionalRecOptions\")" &> rec.log

exitcode=$?

echo "*! Exit code of recCPass0.C: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "recCPass0.C exited with code $exitcode" > validation_error.message
    exit 10
fi

mv syswatch.log syswatch_rec.log

if [ "$2" == "OCDB" ]; then
    echo "*! Reconstruction ran in fake mode to create OCDB.root, exiting quickly now"

    if [ -f OCDB.root ]; then
        echo "* OCDB.root was indeed produced"
    else
        echo "! Error: OCDB.root was NOT generated !!!"
        echo "OCDB.root was not generated" > validation_error.message
        exit 1
    fi
    
    exit 0
fi

echo "* Running AliRoot to make calibration..."
echo executing aliroot -l -b -q -x "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")"
time aliroot -l -b -q -x "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")" &>> calib.log
exitcode=$?

echo "*! Exit code of runCalibTrain.C: $exitcode"

if [ -f ResidualHistos.root ]; then
    mv ResidualHistos.root ResidualTrees.root
fi

echo "*  Running filtering task *"
filtMacro=$ALICE_PHYSICS/PWGPP/macros/runFilteringTask.C
if [ -f $filtMacro ]; then
    echo AliESDs.root > esd.list
    aliroot -l -b -q "${filtMacro}(\"esd.list\",1000,100,\"${ocdbPath}\")" &> filtering.log
else
    echo "no ${filtMacro} ..."
fi

if [ $exitcode -ne 0 ]; then
    echo "runCalibTrain.C exited with code $exitcode" > validation_error.message
    exit 40
fi

exit 0
