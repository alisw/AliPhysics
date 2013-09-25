#!/bin/bash

# Script to run:
#    1. reconstruction
#    2. calibration
#
# Files assumed to be in working directory:
# recCPass1.C          - reconstruction macro
# runCalibTrain.C     - calibration/filtering macro
# Arguments (run locally):
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 

# example:
# runCPass1.sh raw.root  50  104892

#ALIEN setting
# $1 = raw input filename
runNum=`echo $1 | cut -d "/" -f 6 | sed 's/^0*//'`

if [ $# -eq 1 ]; then
    # alien Setup
    nEvents=99999999
    ocdbPath="raw://"
    triggerOptions="?Trigger=kCalibBarrel"
fi

if [ $# -ge 4 ]; then
    # local setup
    nEvents=$2
    runNum=$3
    ocdbPath=$4
    triggerOptions="?Trigger=kCalibBarrel"
fi

if [ $# -eq 5 ]; then
    # local setup in case we specify the trigger mask
    triggerOptions=$5
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
    # several chunks accessed remotely
    CHUNKNAME="collection://wn.xml"
fi

echo "* ************************"
echo "* runCPass1.sh $*"
echo "* Chunk name: $CHUNKNAME"
echo "* Run number: $runNum"
echo "* nEvents: $nEvents"
echo "* runNum: $runNum"
echo "* ocdbPath: $ocdbPath"
echo "* triggerOptions: $triggerOptions"
echo "* ************************"

mkdir Barrel OuterDet

if [ -f Run0_999999999_v3_s0.root ]; then
    echo "* TPC correction file found"

    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/

    for DIR in Barrel OuterDet; do
        mkdir -p $DIR/TPC/Calib/Correction
        ln -s ../../../../TPC/Calib/Correction/Run0_999999999_v3_s0.root $DIR/TPC/Calib/Correction/Run0_999999999_v3_s0.root
    done
fi

echo "* PATH: $PATH"
echo "* LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo

if [ "$2" == "OCDB" ]; then
    export OCDB_SNAPSHOT_CREATE="kTRUE"
    export OCDB_SNAPSHOT_FILENAME="OCDB.root"

    echo "* Running AliRoot to generate the OCDB based on $CHUNKNAME"

    echo "OCDB/recCPass1.C" >&2
    time aliroot -l -b -q recCPass1.C\(\"$CHUNKNAME\"\) &> rec.log
    exitcode=$?
    echo "Exit code: $exitcode"

    echo "* Reconstruction ran in fake mode to create OCDB.root, exiting quickly now"
    touch OCDB.generating.job

    if [ -f OCDB.root ]; then
        echo "* OCDB.root was indeed produced"
    else
        echo "* Error: OCDB.root was NOT generated !!!"
        exit 1
    fi

    exit $exitcode
fi

for COMMON_FILE in wn.xml localOCDBaccessConfig.C AddTaskTPCCalib.C OCDB.root QAtrain_duo.C; do
    if [ -f $COMMON_FILE ]; then
        ln -s ../$COMMON_FILE Barrel/$COMMON_FILE
        ln -s ../$COMMON_FILE OuterDet/$COMMON_FILE
    fi
done

for BARREL_FILE in recCPass1.C runCalibTrain.C; do
    ln -s ../$BARREL_FILE Barrel/$BARREL_FILE
done

for OUTER_FILE in recCPass1_OuterDet.C; do
    ln -s ../$OUTER_FILE OuterDet/$OUTER_FILE
done

####################################   Barrel   #######################################

cd Barrel

echo "* Running AliRoot to reconstruct barrel of $CHUNKNAME"

echo executing aliroot -l -b -q "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")"
time aliroot -l -b -q "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")" &> ../rec.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    exit $exitcode
fi

mv AliESDs.root ../AliESDs_Barrel.root
mv AliESDfriends.root ../AliESDfriends_Barrel.root
mv syswatch.log ../syswatch_rec_Barrel.log

echo "* Running AliRoot to make calibration..."

echo executing aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")"
time aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")" &> ../calib.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    exit $exitcode
fi

mv AliESDfriends_v1.root ../
mv syswatch.log ../syswatch_calib.log

if [ -f QAtrain_duo.C ]; then
    echo "* Running the QA train (barrel) ..."

    echo executing aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"wn.xml\",0,\"$ocdbPath\")"
    time aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"wn.xml\",0,\"$ocdbPath\")" &> ../qa_barrel.log

    exitcode=$?
    echo "Exit code: $exitcode"

    if [ $exitcode -ne 0 ]; then
        exit $exitcode
    fi

    for file in *.stat; do
        mv $file ../$file.qa_barrel
    done
fi

for file in QAresults_barrel.root EventStat_temp_barrel.root AODtpITS.root; do
    if [ -f "$file" ]; then
        mv "$file" ../
    fi
done

####################################   Outer   #######################################

cd ../OuterDet

echo "* Running AliRoot to reconstruct outer of $CHUNKNAME"

echo executing aliroot -l -b -q "recCPass1_OuterDet.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\")"
time aliroot -l -b -q "recCPass1_OuterDet.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\")" &> ../rec_Outer.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    exit $exitcode
fi

mv syswatch.log ../syswatch_rec_Outer.log

if [ -f QAtrain_duo.C ]; then
    echo "* Running the QA train (outer) ..."

    echo executing aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"wn.xml\",0,\"$ocdbPath\")"
    time aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"wn.xml\",0,\"$ocdbPath\")" &> ../qa_outer.log

    exitcode=$?
    echo "Exit code: $exitcode"

    if [ $exitcode -ne 0 ]; then
        exit $exitcode
    fi

    for file in *.stat ; do
        mv $file ../$file.qa_outer
    done
fi

mv AliESDs.root ../AliESDs_Outer.root
mv AliESDfriends.root ../AliESDfriends_Outer.root

for file in QAresults_outer.root EventStat_temp_outer.root; do
    if [ -f "$file" ]; then
        mv "$file" ../
    fi
done

exit $exitcode
