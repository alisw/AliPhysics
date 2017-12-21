#!/bin/bash
# 
# This is the scropt to run the TPC SP calibration procedure as implemented
# in AliTPCDcalibRes
#
# arguments: 
# 1) start time
# 2) end time
# 3) runNumber
# 4) optional his IDs (if missing: all components 0-5)
Usage() {
    echo "Usage: ${0##*/} <minTime> <maxTime> <runNumber> "
    exit
}

if [[ $# -lt 3 ]] ; then Usage ;fi

echo "Arguments: 1 = $1, 2 = $2, 3 = $3"

source $ALICE_ROOT/libexec/scripts/alilog4bash.sh
#

residualFileList="residualFilesList.log" # file containing the whole list of ResidualTrees.root files to be processed
export onGrid=`grep "alien://" $residualFileList | wc -l`

# check if the text file with residual trees is in the working directory
if [[ ! -e $residualFileList ]] ; then
    alilog_info "Error: list of residual trees $residualFileList is not in working directory"
    exit
fi

# check if the drift velocity calibration is in the working directory
if [[ ! -e fitDrift.root ]] ; then
    alilog_info "Error: drift velocity calibration is not in working directory"
    exit
fi

#
export mapStartTime=$1
export mapStopTime=$2
export runNumber=$3

#use of TOF BC or not
export useTOFBC
useTOFBC=${ALIEN_JDL_USETOFBC-$useTOFBC}
#
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/CreateResCalib.C"

#
# we need macro locally to compile it
locMacro=$(basename "$macroName")
if [[ ! -f "$locMacro" ]] ; then cp $macroName ./ ; fi
#
run=$(echo "$runNumber" | sed 's/^0*//')
#
alilog_info "BEGIN Processing for time bin $mapStartTime : $mapStopTime in run $runNumber"
#time aliroot -b -q  $inclMacro ${locMacro}+g\($run,$mapStartTime,$mapStopTime,\"$residualFileList\"\) >& out_${mapStartTime}_${mapStopTime}.log
time aliroot -b -q  $inclMacro ${locMacro}\($run,$mapStartTime,$mapStopTime,\"$residualFileList\"\) >& out_${mapStartTime}_${mapStopTime}.log
alilog_info "END: Processing"
if [[ -f syswatch.log ]] ; then mv syswatch.log syswatch_${mapStartTime}_${mapStopTime}.log ; fi
rm tmpDeltaSect*.root




