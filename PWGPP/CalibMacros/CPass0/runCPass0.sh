#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration 
#
# Files assumed to be in working directory:
# recCPass0.C          - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments (for running on the grid, as called from JDL in central productions):
#    1 - raw data file
#    2 - "OCDB" for creating the OCDB snapshot

# Arguments (local mode - triggered when $# >= 4):
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 
#    4  - OCDBPath
#    5  - optional trigger mask
# example:
# runCPass0.sh raw.root  50  104892 raw://

#ALIEN setting
# $1 = raw input filename
runNum=`echo $1 | cut -d "/" -f 6 | sed 's/^0*//'`
if [ $# -lt 3 ] ; then
  # alien Setup
  nEvents=99999999
  fileName="alien://"$1
  ocdbPath="raw://"
  triggerAlias="?Trigger=kCalibBarrel"
fi;
if [ $# -ge 4 ] ; then
  # local setup
  fileName=$1
  nEvents=$2
  runNum=$3
  ocdbPath=$4
  triggerAlias="?Trigger=kCalibBarrel"
fi
if [ $# -eq 5 ] ; then
  # local setup in case we provide the trigger mask
  # the trigger mask is first stripped of quotation characters
  triggerAlias=${5//\"/}
fi

echo xxxxxxxxxxxxxxxxxxxxxxxxxxx
echo runCPass0.sh Input arguments
echo fileName=$fileName
echo nEvents=$nEvents
echo runNum=$runNum
echo ocdbPath=$ocdbPath
echo triggerAlias=$triggerAlias
echo xxxxxxxxxxxxxxxxxxxxxxxxxxx

if [ -f Run0_999999999_v3_s0.root ]; then
    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/
fi



echo File to be  processed $1
echo Number of events to be processed $nEvents

echo ">>>>>>>>> PATH is..."
echo $PATH
echo ">>>>>>>>> LD_LIBRARY_PATH is..."
echo $LD_LIBRARY_PATH
echo ">>>>>>>>> recCPass0.C is..."
#cat recCPass0.C
echo

echo ">>>>>>> Running AliRoot to reconstruct $1. Run number is $runNum..."

if [ "$2" == "OCDB" ]; then
    echo "Generating OCDB.root only"
    export OCDB_SNAPSHOT_CREATE="kTRUE"
    export OCDB_SNAPSHOT_FILENAME="OCDB.root"
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
echo ""

echo aliroot -l -b -q "recCPass0.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerAlias\")"
time aliroot -l -b -q "recCPass0.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerAlias\")" &> rec.log

exitcode=$?

echo "*! Exit code of recCPass0.C(\"$CHUNKNAME\"): $exitcode"

mv syswatch.log syswatch_rec.log
echo "directory contents:"
ls

if [ "$2" == "OCDB" ]; then
    echo "*! Reconstruction ran in fake mode to create OCDB.root, exiting quickly now"
    touch OCDB.generating.job

    if [ -f OCDB.root ]; then
        echo "* OCDB.root was indeed produced"
    else
        echo "! Error: OCDB.root was NOT generated !!!"
        exit 1
    fi
    exit 0
fi

echo "* Running AliRoot to make calibration..."
echo time aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")"
time aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")" &> calib.log
exitcode=$?

echo "*! Exit code of runCalibTrain.C(\"$runNum\"): $exitcode"

mv syswatch.log syswatch_calib.log

echo ">>>>>>> Extracting system information..."
echo aliroot -b -q $ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeSyswatchCPass0.C\(\"AliESDfriends_v1.root\"\)
aliroot -b -q $ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeSyswatchCPass0.C\(\"AliESDfriends_v1.root\"\)
