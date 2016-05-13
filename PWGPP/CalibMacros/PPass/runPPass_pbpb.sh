#!/bin/bash

echo ALICE_ROOT = $ALICE_ROOT
echo AliROOT = $AliROOT
cp $ALICE_ROOT/.rootrc ~/.rootrc
cp $ALICE_ROOT/.rootrc $HOME
#cat $HOME/.rootrc
export GRID_TOKEN=OK
export XRD_TRANSACTIONTIMEOUT=300

export PRODUCTION_METADATA="$ALIEN_JDL_LPMMETADATA"

echo "* PATH: $PATH"
echo "* LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "* rec.C:"
cat rec.C
echo

timeUsed=0

# run TPC clusterization in separate process before reconstruction
export preclusterizeTPC='1'

ls -l

#ALIEN setting
# $1 = raw input filename
CHUNKNAME="$1"
shift

# second argument could be OCDB
if [ "$1" == "OCDB" ]; then
    echo "* Generating OCDB.root only"
    export OCDB_SNAPSHOT_CREATE="kTRUE"
    export OCDB_SNAPSHOT_FILENAME="OCDB.root"

    shift
fi

# second argument could be SPLIT, then 
# then:
#   $2 - nEvents to reconstruct
#   $3 - runNumber
#   $4 - ocdbPath
if [ "$1" == "SPLIT" ]; then
    nEvents=${2-"-1"}
    runNumber=${3-""}
    ocdbPath=${4-""}
    shift 2
fi

[[ -z $runNumber ]] && runNumber=$(echo $CHUNKNAME | cut -d "/" -f 6 | grep -e '000[0-9][0-9]*' | sed 's/^0*//')
[[ -z $runNumber ]] && runNumber="${CHUNKNAME:2:9}"

echo "runNumber=$runNumber"
echo "ocdbPath=$ocdbPath"
echo "nEvents=$nEvents"
echo "additionalRecOptions=$additionalRecOptions"

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

# Extraction of TPC clusters
if [ -n "$preclusterizeTPC" ] && [ "$OCDB_SNAPSHOT_CREATE" != "kTRUE" ]; then
    echo executing aliroot -l -b -q -x "raw2clust.C(\"$CHUNKNAME\",$nEvents,\"$ocdbPath\")"
    timeStart=`date +%s`
    time aliroot -l -b -q -x "raw2clust.C(\"$CHUNKNAME\",$nEvents,\"$ocdbPath\")" &> clust.log
    exitcode=$?
    timeEnd=`date +%s`
    timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
    delta=$(( $timeEnd-$timeStart ))
    echo "raw2clust:   delta = $delta, timeUsed so far = $timeUsed"
    echo "*! Exit code of raw2clust.C: $exitcode"
    if [ $exitcode -ne 0 ]; then
	echo "raw2clust.C exited with code $exitcode" > validation_error.message
	exit 10
    fi
# Cleanup after the extraction of TPC clusters
    rm -f galice.root AliESDs.root Run*.root Trigger.root QA.root
fi


echo "rec.C" >&2
echo executing aliroot -l -b -q -x "rec.C(\"$CHUNKNAME\",$nEvents,\"$ocdbPath\",\"$additionalRecOptions\")"
timeStart=`date +%s`
time aliroot -l -b -q -x "rec.C(\"$CHUNKNAME\",$nEvents,\"$ocdbPath\",\"$additionalRecOptions\")" &> rec.log
exitcode=$?
timeEnd=`date +%s`
timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
delta=$(( $timeEnd-$timeStart ))
echo "rec:   delta = $delta, timeUsed so far = $timeUsed"

echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "rec.C exited with code $exitcode" > validation_error.message
    exit 10
fi

if [ "$OCDB_SNAPSHOT_CREATE" == "kTRUE" ]; then
    echo "* Reconstruction ran in fake mode to create OCDB.root, exiting quickly now"
    touch OCDB.generating.job

    if [ -f OCDB.root ]; then
        echo "*  OCDB.root was indeed produced"
    else
        echo "!!! Error: OCDB.root was NOT generated !!!"
        echo "OCDB.root was not generated" > validation_error.message
        exit 1
    fi
    
    exit 0
fi

if [ -n "$preclusterizeTPC" ] && [ -f TPC.RecPoints.root ]; then
# clean TPC recpoints
    echo "Removing preproduced TPC recpoints"
    rm TPC.RecPoints.root
fi


echo "* Running AliRoot to generate Tags..."
echo "tag.C" >&2
timeStart=`date +%s`
time aliroot -l -b -q -x tag.C\(\) &> tag.log
exitcode=$?
timeEnd=`date +%s`
timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
delta=$(( $timeEnd-$timeStart ))
echo "tag:   delta = $delta, timeUsed so far = $timeUsed"

echo "Exit code: $exitcode"

if [ $exitcode -ne 0 ]; then
    echo "tag.C exited with code $exitcode" > validation_error.message
    exit 50
fi

for file in *.stat; do
    mv $file $file.rec
done

if [ -f QAtrain.C ]; then
    echo "* Running the QA train..."
    timeStart=`date +%s`
    time aliroot -b -q -x QAtrain.C\($runNumber\) &> qa.log
    exitcode=$?
    timeEnd=`date +%s`
    timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
    delta=$(( $timeEnd-$timeStart ))
    echo "QA:   delta = $delta, timeUsed so far = $timeUsed"
    echo "QAresults.root" >> validation_extrafiles.list

    echo "Exit code: $exitcode"
    
    if [ $exitcode -ne 0 ]; then
        echo "QAtrain.C exited with code $exitcode" > validation_error.message
        exit 100
    fi

    for file in *.stat; do
        mv $file $file.qa
    done
fi

if [ -f QAtrain_duo.C ]; then
    echo "QAresults.root" >> validation_extrafiles.list

    for grp in 0 1 2 3 4; do
        export QAGROUP=$grp
	    echo "* Running QA for tasks group $QAGROUP"
	    echo executing aliroot -b -q -x "QAtrain_duo.C(\"_grp$grp\",$runNumber,0,0,\"$ocdbPath\")"
        echo "QAtrain_duo.C / grp$grp" >&2
	timeStart=`date +%s`
	time aliroot -b -q -x "QAtrain_duo.C(\"_grp$grp\",$runNumber,0,0,\"$ocdbPath\")" &> qa_grp$grp.log
    	exitcode=$?
	timeEnd=`date +%s`
	timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
	delta=$(( $timeEnd-$timeStart ))
	echo "QA $grp:   delta = $delta, timeUsed so far = $timeUsed"
	    echo "Exit code: $exitcode"
	    if [ $exitcode -ne 0 ]; then
            echo "!!! QAtrain_duo.C / barrel group $grp exited with code $exitcode" > validation_error.message
            # put the partial results in the main folder so they are registered, for QA debugging purposes
            exit 100
    	fi
        
        for file in *.stat; do
            mv $file $file.qa_grp$grp
        done
    done
    
    if [ -f mergeQAgroups.C ]; then
        lstQA="lstQA.txt"
	    ls QAresults_grp*.root > $lstQA
	    echo executing aliroot -b -q -x "mergeQAgroups.C(\"$lstQA\",\"QAresults.root\")"
        echo "mergeQAgroups.C" >&2
	timeStart=`date +%s`
	time aliroot -b -q -x "mergeQAgroups.C(\"$lstQA\",\"QAresults.root\")"  &> mergeQA.log
	exitcode=$?
	timeEnd=`date +%s`
	timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
	delta=$(( $timeEnd-$timeStart ))
	echo "QA merge:   delta = $delta, timeUsed so far = $timeUsed"
    else
	    echo "!!! No mergeQAgroups.C macro in current directory, cannot merge..."
    fi
fi

if [ -f AODtrain.C ]; then
    rm -f outputs_valid &>/dev/null
    echo "AliAOD.root" >> validation_extrafiles.list

    echo "* Running the FILTERING train..."
    echo "AODtrain.C" >&2
    timeStart=`date +%s`

    #third argument here is 1 for pbpb
    time aliroot -b -q  -x "AODtrain.C(0,\"$ocdbPath\",1)" &> aod.log
    
    exitcode=$?
    timeEnd=`date +%s`
    timeUsed=$(( $timeUsed+$timeEnd-$timeStart ))
    delta=$(( $timeEnd-$timeStart ))
    echo "aod:   delta = $delta, timeUsed so far = $timeUsed"

    echo "The time spent in aliroot sessions is $timeUsed"
    echo "Exit code: $exitcode"
    
    if [ $exitcode -ne 0 ]; then
        echo "AODtrain.C exited with code $exitcode" > validation_error.message
        exit 200
    fi

    for file in *.stat; do
        mv $file $file.aod
    done
fi

exit 0


