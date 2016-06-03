#!/bin/bash

Usage() {
    echo "Usage: ${0##*/} <xml collection of residual trees> <runNum1> [runNum2] ..."
    echo "where runNum's are the runs to be processed; at least one MUST be provided"
    exit
}

##############################################################################
extractFileNamesFromXMLCollection()
{
    egrep turl|sed 's|^.*turl\s*=\s*"\s*\([a-zA-Z]*://.*\.root\).*$|\1|g'
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
[[ $# -ne 2 ]] &&  Usage && exit

echo Starting...
residualFileList=residualFilesList.log
runList=runList.log

if [[ $1 == *.xml ]]; then  # case of xml collection: we need the list of filenames
#echo "the XML file to be parsed is:"
#cat $1
    cat $1 | extractFileNamesFromXMLCollection > $residualFileList
    echo
else
    if [[ $1 != $residualFileList ]]; then  # The input file list is not in a file with the filename that we want. We need this for the following step (processTimeBin)
	ln -s $1 $residualFileList
    fi
fi

shift 1
[[ -f $runList ]] && rm $runList

nruns=$#
for runNum in "$@"
do
    echo $runNum >> $runList
done


echo "The runs to be processed are"
echo `cat $runList`

export onGrid=`grep "alien://" $residualFileList | wc -l`
echo onGrid=$onGrid
export treeCacheSize=100000000
export autoCacheSize=0
echo "onGrid=$onGrid, treeCacheSize=$treeCacheSize, autoCacheSize=$autoCacheSize"

echo sourcing alilog4bash.sh
source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh  

loadLibMacro="$ALICE_PHYSICS/PWGPP/CalibMacros/CPass1/LoadLibraries.C"
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procResidData.C"

extractEnvVars

curdir=`pwd`
#
# we need macro locally to compile it
locMacro=$(basename "$macroName")
[[ ! -f "$locMacro" ]] && cp $macroName ./ 
[[ ! -f "$locMacro" ]] && echo "did not find $locMacro" && exit 1

residFilesRun="residual.list"

# if >1 run is requested, to avoid ovewriting, results for each run will appear in separate directory
runDirPref="res_"
if [ "$nruns" -ne 1 ] ; then 
    runDirPref="res_"
    alilog_info "***"
    alilog_info "*** ATTENTION: $nruns submitted for processing"
    alilog_info "*** Output of each run will be written in"
    alilog_info "*** $pref_runNumber directory"
    alilog_info "***"
fi

# 1.)  define time bins and extract vdrift
#
wdir=`pwd` 
for arun in `cat $runList`; do
    alilog_info "BEGIN: Processing run $arun"
    run=`echo $arun| sed s_000__`
    echo "run=$run"
    export runNumber=$run
    runDir=$curdir
    if [ -z "$runDirPref" ] ; then
	runDir="$curdir/$runDirPref_$run"
	[ ! -d /tmp/mydir  ] && mkdir $runDir
	echo "Results will be stored in $runDir directory"
	cd $runDir
    fi
    cat "$curdir/$residualFileList" | egrep $arun > $residFilesRun ;	
    #
    mode=0; 
    alilog_info "BEGIN: 0.) Submit in mode $mode to get time/drift info"
    time aliroot -b -q $inclMacro $loadLibMacro $curdir/${locMacro}+g\("$mode","$run",0,-1,\""$residFilesRun"\"\) >& submitTimeDrift_$arun.log	
    rm tmpDriftTree.root
    alilog_info "END: 0.) Submit int mode $mode to get time/vdrift info"
    #
    mode=1
    alilog_info "BEGIN: 1.) Submit in mode $mode to create OCDB entry for drift velocity"
    time aliroot -b -q $inclMacro $loadLibMacro $curdir/${locMacro}+g\("$mode","$run"\) >& ocdb_vdrift_$arun.log
    alilog_info "END: 1.)  Submit in mode $mode to create OCDB entry for drift velocity"
    alilog_info "END: Processing run $arun"
    #
    timeBinsFile="timeBins.log"
    alilog_info "BEGIN: 2.) Writing time bins info to $timeBinsFile"
    # for mac: egrep instead of grep, and awk instead of gawk
    gminTime=`cat submitTimeDrift_$arun.log  | egrep StatInfo.minTime | awk '{print $2}' | tail -n 1`
    gmaxTime=`cat submitTimeDrift_$arun.log  | egrep StatInfo.maxTime | awk '{print $2}' | tail -n 1`
    gNTimeBins=`cat submitTimeDrift_$arun.log  | egrep StatInfo.NTBin | awk '{print $2}' | tail -n 1`
    gTimeDelta=`cat submitTimeDrift_$arun.log  | egrep StatInfo.TBinDuration | awk '{print $2}' | tail -n 1`
    tstart=$gminTime
    tend=$gminTime
    alilog_info "Time bins for run $run: $gNTimeBins of $gTimeDelta seconds"
    for ((itime=1;itime<=$gNTimeBins;itime+=1));  do   # time loop
	tend=`echo $gmaxTime $gminTime $gNTimeBins $itime | awk '{printf "%d\n", $2+($1-$2)/$3*$4}'`
	alilog_info "BEGIN: Writing time bin $itime for run $arun : $tstart - $tend"
	echo $tstart $tend $arun >> $timeBinsFile
	alilog_info "END: Writing time bin $itime for run $arun : $tstart - $tend"
	tstart=$tend
    done;
    #   
    alilog_info "End: 2.) Writing time bins info"
    cd $curdir
    alilog_info "End: Processing run $arun"
done

alilog_info "END: VDrift and time bins definition"
