#!/bin/bash

# script will extract map for already preprocessed data, starting from
# alitpcdcalibres.root in current directory
# arguments: 
# 1) start time
# 2) end time
# 3) run number
# 4) optional number of tracks for closure test (if requested)


##############################################################################
Usage() {
    echo "Usage: ${0##*/} <minTime> <maxTime> <runNumber> [ntracks_closure_test]"
    exit
}

##############################################################################
extractEnvVars()
{
    # variables to calculate vdrfit
    export driftDeltaT=${ALIEN_JDL_DRIFTDELTAT-$driftDeltaT}
    export driftSigmaT=${ALIEN_JDL_DRIFTSIGMAAT-$driftSigmaT}
    
    export targetOCDBDir=${ALIEN_JDL_TARGETOCDBDIR-$targetOCDBDir}
    if [ -z "$targetOCDBDir" ] ; then targetOCDBDir="local://`pwd`" ; fi 
    
    export distTimeBin=${ALIEN_JDL_DISTTIMEBIN-$distTimeBin}  # nominal time bin 
    
    export distTimeBinPrec=${ALIEN_JDL_DISTTIMEBINPREC-$distTimeBinPrec}; # nominal time bin relative precision (can play within this factor)
    
    # do we use TOF BC validataion
    export useTOFBC=${ALIEN_JDL_USETOFBC-$useTOFBC}
    
    # which detectors to use
    export distUseDet=${ALIEN_JDL_DISTUSEDET-$distUseDet}
    
    export distNBinsZ=${ALIEN_JDL_DISTNBINSZ-$distNBinsZ}
    export distNBinsY=${ALIEN_JDL_DISTNBINSY-$distNBinsY}
    
    export distMaxTracks=${ALIEN_JDL_DISTMAXTRACKS-$distMaxTracks}
    export distMinTracks=${ALIEN_JDL_DISTMINTRACKS-$distMinTracks}

    export distMinValidVoxPerRow=${ALIEN_JDL_DISTMINVALIDVOXPERROW-$distMinValidVoxPerRow}
    
    export distKernelType=${ALIEN_JDL_DISTKERNELTYPE-$distKernelType}
    export distKernelWX=${ALIEN_JDL_DISTKERNELWX-$distKernelWX}
    export distKernelWY=${ALIEN_JDL_DISTKERNELWY-$distKernelWY}
    export distKernelWZ=${ALIEN_JDL_DISTKERNELWZ-$distKernelWZ}
    
    export distKernelPol2X=${ALIEN_JDL_DISTKERNELPOL2X-$distKernelPol2X}
    export distKernelPol2Y=${ALIEN_JDL_DISTKERNELPOL2Y-$distKernelPol2Y}
    export distKernelPol2Z=${ALIEN_JDL_DISTKERNELPOL2Z-$distKernelPol2Z}
    
    # number of tracks for closure test, hardly will be done of the grid
    export distNTracksClosureTest=${ALIEN_JDL_DISTNTRACKSCLOSURETEST-$distNTracksClosureTest}
    #
    echo ""
    echo "Listing all env.vars"
    printenv
    echo ""
}
##############################################################################


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

extractEnvVars

run=$(echo "$runNumber" | sed 's/^0*//')
alilog_info "BEGIN Processing for time bin $mapStartTime : $mapStopTime in run $runNumber"
mode=2
time aliroot -b -q  $inclMacro $loadLibMacro ${locMacro}+g\($mode,$run,$mapStartTime,$mapStopTime,\"\"\) >& out_${mapStartTime}_${mapStopTime}.log
alilog_info "END: Processing"
[[ -f syswatch.log ]] && mv syswatch.log syswatch_${mapStartTime}_${mapStopTime}.log 
rm tmpDeltaSect*.root

if [ -n "$ntcloseTest" ] ; then
    alilog_info "BEGIN Closure test with $ntcloseTest tracks for time bin $mapStartTime : $mapStopTime in run $runNumber"
    mode=3
    time aliroot -b -q  $inclMacro $loadLibMacro ${locMacro}+g\($mode,$run,$mapStartTime,$mapStopTime,\"\"\) >& closure_${mapStartTime}_${mapStopTime}.log
    alilog_info "END: Processing"
fi
