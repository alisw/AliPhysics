#!/bin/bash

# script will extract map for already preprocessed data, starting from
# alitpcdcalibres.root in current directory
# arguments: 
# 1) start time
# 2) end time
# 3) run number
# 4) optional number of tracks for closure test (if requested)

Usage() {
    echo "Usage: ${0##*/} <minTime> <maxTime> <runNumber> [ntracks_closure_test]"
    exit
}

[[ $# -lt 3 ]] &&  Usage && exit
echo "Arguments: 1 = $1, 2 = $2, 3 = $3"

source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh  

prepobj="alitpcdcalibres.root"

# check if the preprocessed object exists
if [[ ! -e $prepobj ]] ; then
    alilog_info "Error: preprocessed object $prepobj is not in working directory"
    exit
fi

export mapStartTime=$1
export mapStopTime=$2
export runNumber=$3
export distNTracksClosureTest
if [[ $# -gt 3 ]] ; then distNTracksClosureTest=$4 ;fi

loadLibMacro="$ALICE_PHYSICS/PWGPP/CalibMacros/CPass1/LoadLibraries.C"
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procResidData.C"
locMacro=$(basename "$macroName")
[[ ! -f "$locMacro" ]] && cp $macroName ./
[[ ! -f "$locMacro" ]] && echo "did not find $locMacro" && exit 1

run=$(echo "$runNumber" | sed 's/^0*//')
alilog_info "BEGIN Processing for time bin $mapStartTime : $mapStopTime in run $runNumber"
mode=2
time aliroot -b -q  $inclMacro $loadLibMacro ${locMacro}+g\($mode,$run,$mapStartTime,$mapStopTime,\"\"\) >& out_${mapStartTime}_${mapStopTime}.log
alilog_info "END: Processing"
[[ -f syswatch.log ]] && mv syswatch.log syswatch_${mapStartTime}_${mapStopTime}.log 
rm tmpDeltaSect*.root

if [ -z "$ntcloseTest" ] ; then
    alilog_info "BEGIN Closure test with $ntcloseTest tracks for time bin $mapStartTime : $mapStopTime in run $runNumber"
    mode=3
    time aliroot -b -q  $inclMacro $loadLibMacro ${locMacro}+g\($mode,$run,$mapStartTime,$mapStopTime,\"\"\) >& closure_${mapStartTime}_${mapStopTime}.log
    alilog_info "END: Processing"
fi
