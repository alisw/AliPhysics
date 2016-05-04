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
# runCPass1_pbpb2015_pass2.sh raw.root  50  104892

#ALIEN setting
# $1 = raw input filename
tmpName=$(basename "$1")
runNumF="${tmpName:2:9}"
runNum=`echo "$runNumF" |  sed 's/^0*//'`

#default alias for CPass1
defAlias="kCalibBarrelMB"

#optionallly skip Outer pass
export skipOuter='1'

# Exporting variable to define that we are in CPass1 to be used in reconstruction
export CPass='1'

export PRODUCTION_METADATA="$ALIEN_JDL_LPMMETADATA"

# Set memory limits to a value lower than the hard site limits to at least get the logs of failing jobs
ulimit -S -v 3500000

# run TPC clusterization in separate process before reconstruction
export preclusterizeTPC='1'

if [ $# -eq 1 ]; then
    # alien Setup
    nEvents=99999999
    ocdbPath="raw://"
    # use the value passed by LPM, or by default use the kCalibBarrel alias
    triggerOptions=${ALIEN_JDL_TRIGGERALIAS-?Trigger=$defAlias}
    #triggerOptions="?Trigger=kPhysicsAll"
fi

if [ $# -ge 4 ]; then
    # local setup
    nEvents=$2
    runNum=$3
    ocdbPath=$4
    triggerOptions="?Trigger=$defAlias"
fi

if [ $# -eq 5 ]; then
    # local setup in case we specify the trigger mask
    triggerOptions=$5
fi

# ===| TPC default values |===================================================
# can be overwritten by config file, or JDL below
# JDL will overwrite the config file
#
# ---| gain calibration |-----------------------------------------------------
# default in CPass1 is combined calibration + residual QA,
#    for number convention see AliTPCPreprocessorOffline::EGainCalibType
#
export TPC_CPass1_GainCalibType=3

# ===| TPC JDL overwrites |===================================================
#
export TPC_CPass1_GainCalibType=${ALIEN_JDL_TPC_CPASS1_GAINCALIBTYPE-$TPC_CPass1_GainCalibType}

echo "TPC_CPass1_GainCalibType=${TPC_CPass1_GainCalibType}" | tee -a calib.log

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
if [ -n  "$skipOuter" ]; then
    echo "Outer pass will be skept"
fi

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
    touch OCDB.generating.job

    echo "* Running AliRoot to generate the OCDB based on $CHUNKNAME"

    echo "OCDB/recCPass1.C" >&2
    time aliroot -l -b -q -x recCPass1.C\(\"$CHUNKNAME\"\) &> rec.log
    exitcode=$?
    echo "Exit code: $exitcode"
    
    if [ $exitcode -ne 0 ]; then
        echo "recCPass1.C for OCDB snapshot exited with code $exitcode" > validation_error.message
        exit 10
    fi

    echo "* Reconstruction ran in fake mode to create OCDB.root, exiting quickly now"

    if [ -f OCDB.root ]; then
        echo "* OCDB.root was indeed produced"
    else
        echo "* Error: OCDB.root was NOT generated !!!"
        echo "OCDB.root was not generated" > validation_error.message
        exit 1
    fi

    exit 0
fi

for COMMON_FILE in wn.xml localOCDBaccessConfig.C AddTaskTPCCalib.C AddTaskTRDCalib.C OCDB.root QAtrain_duo.C mergeQAgroups.C; do
    if [ -f $COMMON_FILE ]; then
        ln -s ../$COMMON_FILE Barrel/$COMMON_FILE
        ln -s ../$COMMON_FILE OuterDet/$COMMON_FILE
    fi
done

for BARREL_FILE in recCPass1.C runCalibTrain.C raw2clust.C; do
    ln -s ../$BARREL_FILE Barrel/$BARREL_FILE
done

for OUTER_FILE in recCPass1_OuterDet.C; do
    ln -s ../$OUTER_FILE OuterDet/$OUTER_FILE
done

####################################   Barrel   #######################################

cd Barrel

echo "* Running AliRoot to reconstruct barrel of $CHUNKNAME"

# Extraction of TPC clusters
if [ -n "$preclusterizeTPC" ]; then
    echo executing aliroot -l -b -q -x "raw2clust.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")"
    time aliroot -l -b -q -x "raw2clust.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")" &> clust.log
    exitcode=$?
    echo "*! Exit code of raw2clust.C: $exitcode"
    if [ $exitcode -ne 0 ]; then
	echo "raw2clust.C exited with code $exitcode" > validation_error.message
	exit 10
    fi
# Cleanup after the extraction of TPC clusters
    rm -f galice.root AliESDs.root Run*.root Trigger.root QA.root
fi

echo executing aliroot -l -b -q -x "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")"
time aliroot -l -b -q -x "recCPass1.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\", \"$triggerOptions\")" &> ../rec.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "recCPass1.C exited with code $exitcode" > ../validation_error.message
    exit 10
fi

if [ -n "$preclusterizeTPC" ] && [ -f TPC.RecPoints.root ]; then
# clean TPC recpoints
    echo "Removing preproduced TPC recpoints"
    rm TPC.RecPoints.root
    mv clust.log ../clust.log
fi

mv syswatch.log ../syswatch_rec_Barrel.log

echo "* Running AliRoot to make calibration..."

echo executing aliroot -l -b -q -x "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")"
time aliroot -l -b -q -x "runCalibTrain.C($runNum,\"AliESDs.root\",\"$ocdbPath\")" &>> ../calib.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "runCalibTrain.C exited with code $exitcode" > ../validation_error.message
    exit 40
fi

mv syswatch.log ../syswatch_calib.log

if [ -f ResidualHistos.root ]; then
    mv ResidualHistos.root ../ResidualTrees.root
fi
 
echo "*  Running filtering task for barrel *"
filtMacro=$ALICE_PHYSICS/PWGPP/macros/runFilteringTask.C
if [ -f $filtMacro ]; then
    echo AliESDs.root > esd.list
    aliroot -l -b -q "${filtMacro}(\"esd.list\",10000,1000,\"${ocdbPath}\")" &> filtering.log
else
    echo "no ${filtMacro} ..."
fi
if [ -f filtering.log ]; then
    mv filtering.log ../filtering.log
fi
#
if [ -f QAtrain_duo.C ]; then
    echo "* Running the QA train (barrel) ..."

#    echo executing aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"$ocdbPath\")"
#    time aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNum,\"$ocdbPath\")" &> ../qa_barrel.log

    for grp in 0 1 2 3 4 
    do
	export QAGROUP=$grp
	echo running QA for tasks group $QAGROUP
	echo executing aliroot -b -q -x "QAtrain_duo.C(\"_barrel_grp$grp\",$runNum,0,0,\"$ocdbPath\")"
	time aliroot -b -q -x "QAtrain_duo.C(\"_barrel_grp$grp\",$runNum,0,0,\"$ocdbPath\")" &> ../qa_barrel_grp$grp.log
#
	exitcode=$?
	echo "Exit code: $exitcode"
	if [ $exitcode -ne 0 ]; then
            echo "QAtrain_duo.C / barrel group $grp exited with code $exitcode" > ../validation_error.message
        # put the partial results in the main folder so they are registered, for QA debugging purposes
            mv AliESDs.root ../AliESDs_Barrel.root
            mv AliESDfriends.root ../AliESDfriends_Barrel.root
            exit 100
	fi
#
    done
#
    if [ -f mergeQAgroups.C ]; then
	lstQAbarrel="lstQAbarrel.txt"
	ls QAresults_barrel_grp*.root > $lstQAbarrel
	echo executing aliroot -b -q -x "mergeQAgroups.C(\"$lstQAbarrel\",\"QAresults_barrel.root\")"
	time aliroot -b -q -x "mergeQAgroups.C(\"$lstQAbarrel\",\"QAresults_barrel.root\")"  &> ../mergeQA_barrel.log
    else
	echo "no mergeQAgroups.C macro in current directory, cannot merge..."
    fi
 #  

    for file in *.stat; do
        mv $file ../$file.qa_barrel
    done
fi

mv AliESDs.root ../AliESDs_Barrel.root
mv AliESDfriends.root ../AliESDfriends_Barrel.root

for file in FilterEvents_Trees.root AliESDfriends_v1.root QAresults_barrel.root EventStat_temp_barrel_grp*.root AODtpITS.root Run*.Event*_*.ESD.tag.root TOFcalibTree.root T0AnalysisTree.root CalibObjects.root; do
    if [ -f "$file" ]; then
        mv "$file" ../
    fi
done

# cleanup of barrel
cd ../
mv EventStat_temp_barrel_grp0.root EventStat_temp_barrel.root
rm EventStat_temp_*_grp*.root

if [ -n  "$skipOuter" ]; then
  exit 0
fi

####################################   Outer   #######################################

cd OuterDet

echo "* Running AliRoot to reconstruct outer of $CHUNKNAME"

echo executing aliroot -l -b -q -x "recCPass1_OuterDet.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\")"
time aliroot -l -b -q -x "recCPass1_OuterDet.C(\"$CHUNKNAME\", $nEvents, \"$ocdbPath\")" &> ../rec_Outer.log
exitcode=$?
echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "recCPass1_OuterDet.C exited with code $exitcode" > ../validation_error.message
    exit 11
fi

mv syswatch.log ../syswatch_rec_Outer.log

if [ -f QAtrain_duo.C ]; then
    echo "* Running the QA train (outer) ..."

#    echo executing aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"$ocdbPath\")"
#    time aliroot -b -q "QAtrain_duo.C(\"_outer\",$runNum,\"$ocdbPath\")" &> ../qa_outer.log

    for grp in 0 4 
    do
	export QAGROUP=$grp
	echo running QA for tasks group $QAGROUP
	echo executing aliroot -b -q -x "QAtrain_duo.C(\"_outer_grp$grp\",$runNum,0,0,\"$ocdbPath\")"
	time aliroot -b -q -x "QAtrain_duo.C(\"_outer_grp$grp\",$runNum,0,0,\"$ocdbPath\")" &> ../qa_outer_grp$grp.log
#
	exitcode=$?
	echo "Exit code: $exitcode"
	if [ $exitcode -ne 0 ]; then
            echo "QAtrain_duo.C  / outer group $grp exited with code $exitcode" > ../validation_error.message
        # put the partial results in the main folder so they are registered, for QA debugging purposes
            mv AliESDs.root ../AliESDs_Outer.root
            mv AliESDfriends.root ../AliESDfriends_Outer.root
            exit 101
	fi
    done
#
    lstQAouter="lstQAouter.txt"
    ls QAresults_outer_grp*.root > $lstQAouter
    echo executing aliroot -b -q -x "mergeQAgroups.C(\"$lstQAouter\",\"QAresults_outer.root\")"
    time aliroot -b -q -x "mergeQAgroups.C(\"$lstQAouter\",\"QAresults_outer.root\")"  &> ../mergeQA_outer.log
 #   

    for file in *.stat ; do
        mv $file ../$file.qa_outer
    done
fi

mv AliESDs.root ../AliESDs_Outer.root
mv AliESDfriends.root ../AliESDfriends_Outer.root

for file in QAresults_outer.root EventStat_temp_outer_grp*.root; do
    if [ -f "$file" ]; then
        mv "$file" ../
    fi
done

# cleanup of outer
cd ..
mv EventStat_temp_outer_grp0.root EventStat_temp_outer.root
rm EventStat_temp_*_grp*.root

exit 0
