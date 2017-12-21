#!/bin/bash
# 
# Shell script to create a time bin defintion for correction/distortion maps
#   
# $ALICE_PHYSICS/PWGPP/TPC/CalibMacros/defineTimeBins.sh
#
# arguments:  xml collection or residualList.log file

##############################################################################
Usage() {
    echo "Usage: ${0##*/} <xml collection of residual trees> <runNum1> [runNum2] ..."
    echo "where runNum's are the runs to be processed; at least one MUST be provided"
    exit
}

extractFileNamesFromXMLCollection()
{
    egrep turl|sed 's|^.*turl\s*=\s*"\s*\([a-zA-Z]*://.*\.root\).*$|\1|g'
}

##############################################################################

if [[ $# -ne 2 ]] ; then Usage ;fi

echo Starting...
residualFileList=residualFilesList.log
runList=runList.log

if [[ $1 == *.xml ]]; then  # case of xml collection: we need the list of filenames
#echo "the XML file to be parsed is:"
#cat $1
    cat $1 | extractFileNamesFromXMLCollection > $residualFileList
    echo
else
    if [[ $1 != residualFilesList.log ]]; then  # The input file list is not in a file with the filename that we want. We need this for the following step (processTimeBin)
	ln -s $1 $residualFileList
    fi
fi

shift 1
if [[ -f $runList ]] ; then rm $runList ; fi
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
source $ALICE_ROOT/libexec/alilog4bash.sh

loadLibMacro="$ALICE_PHYSICS/PWGPP/CalibMacros/CPass1/LoadLibraries.C"
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C"
#
# we need macro locally to compile it
locMacro=$(basename "$macroName")
if [[ ! -f "$locMacro" ]] ; then cp $macroName ./ ; fi
# 
residFilesRun="residual.list"

# variables to calculate vdrfit
export driftDeltaT=120;     # make drift update every 120 seconds
export driftSigmaT=600;     # smoothing  sigma for drift calibration

# OCDB directory where to export the drift velocity
export targetOCDBDir="local://`pwd`"
targetOCDBDir=${ALIEN_JDL_TARGETOCDBDIR-$targetOCDBDir}

#use of TOF BC or not
export useTOFBC
useTOFBC=${ALIEN_JDL_USETOFBC-$useTOFBC}
#
# 1.) Submit query  to get  time dependent info
#
wdir=`pwd` 
for arun in `cat $runList`; do
    alilog_info "BEGIN: Processing run $arun"
    run=`echo $arun| sed s_000__`
    echo "run=$run"
    export runNumber=$run
    cat $residualFileList | egrep $arun > $residFilesRun ;	
    alilog_info "BEGIN: 1.) Submit query to get time dependent info"
    time aliroot -b -q $inclMacro $loadLibMacro ${locMacro}+\(5,$run\) >& submitTime_$arun.log	
    alilog_info "END: 1.) Submit query to get time dependent info"
    alilog_info "BEGIN: 1.) Submit query to calibrate drift velocity"
    time aliroot -b -q $inclMacro $loadLibMacro ${locMacro}+\(6,$run\) >& submitDrift_$arun.log
    alilog_info "END: 1.) Submit query to calibrate drift velocity"
    alilog_info "BEGIN: 1.) Create OCDB entry for drift velocity"
    time aliroot -b -q $inclMacro $loadLibMacro ${locMacro}+\(7,$run\) >& ocdb_vdrift_$arun.log
    alilog_info "END: 1.)  Create OCDB entry for drift velocity"
    alilog_info "END: Processing run $arun"
done
#
# 2.) extract bin info
#
alilog_info "BEGIN: Define time bins"
timeDelta=2400 # the width of the time bin will be at minimum 600 s, and at maximum 1200 s 
timeDeltaMinFrac="0.7"  # if possible, try to avoid time-bins smaller than timeDelta*timeDeltaMinFrac


timeDeltaMin=`echo "($timeDelta*$timeDeltaMinFrac)/1" | bc`

timeBinsFile="timeBins.log"

for arun in `cat $runList`; do
#
    if [[ -f $timeBinsFile ]] ; then rm -rf $timeBinsFile ; fi
    gminTime=`cat submitTime_$arun.log  | egrep StatInfo.minTime | awk '{print $2}' | tail -n 1`	 # for mac: egrep instead of grep, and awk instead of gawk
    gmaxTime=`cat submitTime_$arun.log  | egrep StatInfo.maxTime | awk '{print $2}' | tail -n 1`   # for mac: egrep instead of grep, and awk instead of gawk

    # define number of bins with approx timeDelta duration
    nTimeBins=$(( (gmaxTime-gminTime)/timeDelta+1 ))
    timeDeltaRun=$(( (gmaxTime-gminTime)/nTimeBins ))

    if (( $timeDeltaRun < $timeDeltaMin )) ; then # try to avoid too short time intervals
	if (( $nTimeBins > 1 )) ; then
	    ((nTimeBins--)) ;
	    timeDeltaRun=$(( (gmaxTime-gminTime)/nTimeBins )) ;
	fi
    fi

    alilog_info "BEGIN: Processing run $arun"
    alilog_info "Split run $arun ($gminTime : $gmaxTime) to $nTimeBins of $timeDeltaRun duration"

    tstart=$gminTime
    tend=$gminTime
    
    for ((itime=1;itime<=$nTimeBins;itime+=1));  do   # time loop
	tend=`echo $gmaxTime $gminTime $nTimeBins $itime | awk '{printf "%d\n", $2+($1-$2)/$3*$4}'`
	alilog_info "BEGIN: Writing time bin $itime for run $arun : $tstart - $tend"
	echo $tstart $tend $arun >> $timeBinsFile
	alilog_info "END: Writing time bin $itime for run $arun : $tstart - $tend"
	tstart=$tend
    done;
    alilog_info "END:Processing run $arun"
    
done
alilog_info "END: Define time bins"
#



