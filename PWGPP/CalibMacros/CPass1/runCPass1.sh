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
#    4  - OCDBPath
#    5  - optional trigger mask
# example:
# runCPass1.sh raw.root  50  104892 raw://

#ALIEN setting
# $1 = raw input filename
runNum=`echo $1 | cut -d "/" -f 6 | sed 's/^0*//'`
if [ $# -eq 1 ] ; then
  # alien Setup
  nEvents=99999999
  fileName="alien://"$1
  ocdbPath="raw://"
  triggerOptions="?Trigger=kCalibBarrel"
fi;
if [ $# -ge 4 ] ; then
  # local setup
  fileName=$1
  nEvents=$2
  runNum=$3
  ocdbPath=$4
  triggerOptions="?Trigger=kCalibBarrel"
fi
if [ $# -eq 5 ] ; then
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

echo xxxxxxxxxxxxxxxxxxxxxxxxxxx
echo runCPass1.sh Input arguments
echo fileName=$fileName
echo nEvents=$nEvents
echo runNum=$runNum
echo ocdbPath=$ocdbPath
echo triggerOptions=$triggerOptions
echo xxxxxxxxxxxxxxxxxxxxxxxxxxx
echo "* ************************"
echo "* runCPass1.sh $*"
echo "* Chunk name: $CHUNKNAME"
echo "* Run number: $runNum"
echo "* ************************"

if [ -f Run0_999999999_v3_s0.root ]; then
    echo "* TPC correction file found"

    mkdir -p TPC/Calib/Correction
    mv Run0_999999999_v3_s0.root TPC/Calib/Correction/

    mkdir -p Barrel/TPC/Calib/Correction
    ln -s ../../../../TPC/Calib/Correction/Run0_999999999_v3_s0.root Barrel/TPC/Calib/Correction/Run0_999999999_v3_s0.root

    mkdir -p OuterDet/TPC/Calib/Correction
    ln -s ../../../../TPC/Calib/Correction/Run0_999999999_v3_s0.root OuterDet/TPC/Calib/Correction/Run0_999999999_v3_s0.root
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

mkdir Barrel OuterDet

[[ -f localOCDBaccessConfig.C ]] && cp localOCDBaccessConfig.C Barrel
[[ -f localOCDBaccessConfig.C ]] && cp localOCDBaccessConfig.C OuterDet

cp recCPass1.C Barrel/
cp runCalibTrain.C Barrel/
cp QAtrain_duo.C Barrel/

cp recCPass1_OuterDet.C OuterDet/
cp QAtrain_duo.C OuterDet/

if [ -f wn.xml ]; then
    cp wn.xml Barrel/
    cp wn.xml OuterDet/
fi

if [ -f OCDB.root ]; then
    ln -s ../OCDB.root Barrel/OCDB.root
    ln -s ../OCDB.root OuterDet/OCDB.root
fi

####################################   Barrel   #######################################

cd Barrel

echo "* Running AliRoot to reconstruct barrel of $CHUNKNAME"

echo "Barrel/recCPass1.C" >&2
echo executing time aliroot -l -b -q "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")"
time aliroot -l -b -q "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")" &> ../rec_Barrel.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    exit $exitcode
fi

mv syswatch.log ../syswatch_rec_Barrel.log

echo "* Running AliRoot to make calibration..."

echo "Barrel/recCalibTrain.C" >&2
echo executing time aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")"
time aliroot -l -b -q "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")" &> ../calib.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    exit $exitcode
fi

mv syswatch.log ../syswatch_calib.log

if [ -f QAtrain_duo.C ]; then
    echo "* Running the QA train (barrel) ..."

    echo "Barrel/QAtrain_duo.C" >&2
    echo executing time aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"$ocdbPath\")"
    time aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"$ocdbPath\")" &> ../qa_barrel.log
    exitcode=$?
    echo "Exit code: $exitcode"

    if [ $exitcode -ne 0 ]; then
        exit $exitcode
    fi

    for file in *.stat; do
        mv $file ../$file.qa_barrel
    done
fi

mv AliESDs.root ../AliESDs_Barrel.root
mv AliESDfriends.root ../AliESDfriends_Barrel.root

for file in AliESDfriends_v1.root QAresults_barrel.root EventStat_temp_barrel.root AODtpITS.root; do
    if [ -f "$file" ]; then
        mv "$file" ../
    fi
done

####################################   Outer   #######################################

cd ../OuterDet

echo "* Running AliRoot to reconstruct outer of $CHUNKNAME"

echo "OuterDet/recCPass1_OuterDet" >&2
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

    echo "OuterDet/QAtrain_duo.C" >&2
    echo executing time aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"$ocdbPath\")"
    time aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"$ocdbPath\")" &> ../qa_outer.log
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

echo ">>>>>>> Extracting system information..."
echo executing aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/makeSyswatchCPass1.C(\"AliESDfriends_v1.root\")"
aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/makeSyswatchCPass1.C(\"AliESDfriends_v1.root\")"
