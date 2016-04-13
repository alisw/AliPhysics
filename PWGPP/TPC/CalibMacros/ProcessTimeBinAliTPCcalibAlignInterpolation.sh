#!/bin/bash
#
# This is the script to run the TPC SP calibration procedure as implemented
# in AliTPCcalibAlignInterpolation
#
# arguments: 
# 1) start time
# 2) end time
# 3) runNumber
# 4) optional his IDs (if missing: all components 0-5)
Usage() {
    echo "Usage: ${0##*/} <minTime> <maxTime> <runNumber> [hisID1] [hisID2] ... "
    echo "where hisID is components to produce: 0 - 5 (all components by default)"
    exit
}

if [[ $# -lt 3 ]] ; then Usage ;fi

echo "Arguments: 1 = $1, 2 = $2, 3 = $3, 4 = $4, 5 = $5, 6 = $6, 7 = $7"

source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh  
#
secStep=18  # number of sectors to process in one go
hisMin=0
hisMax=5

residualFileList="residualFilesList.log" # file containing the whole list of ResidualTrees.root files to be processed
export onGrid=`grep "alien://" $residualFileList | wc -l`
resList="residual.list" # file containing the list of ResidualTrees.root files to be processed for the current run only

# check if the text file with residual trees is in the working directory
if [[ ! -e $residualFileList ]] ; then
    alilog_info "Error: list of residual trees $residualFileList is not in working directory"
    exit
fi

cat $residualFileList | egrep $3 > $resList 

ln -s fitDrift_$3.root fitDrift.root

minNumberOfEntriesNDFit=10 # check this number
treeNames=( 'deltaRPhiTPCITS' 'deltaRPhiTPCITSTRDDist'  'deltaRPhiTPCITSTOFDist'  'deltaZTPCITSDist'  'deltaZTPCITSTRDDist' 'deltaZTPCITSTOFDist' );
#
# check if the text file with per run residual trees is in the working directory
if [[ ! -e $resList ]] ; then
    alilog_info "Error: list of per run residual trees $resList is not in working directory"
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
export minNumberOfEntriesNDFit
export treeCacheSize=100000000
export autoCacheSize=0
shift 3;
#
# check if the list of map components is provided at the command line
#
nhis=0
for his in "$@"
do
    if [[ $his -ge 0 ]] && [[ $his -le 5 ]] ; then
	hisList[$nhis]=$his
	((nhis+=1))
    else
	alilog_info "Error: map components $his is not limited to $hisMin ... $hisMax"
	exit
    fi
done
#
nhis=${#hisList}
#
if [[ $nhis -lt 1 ]] ; then
    hisList=`seq $hisMin $hisMax`
    hisMin=$3
    hisMax=$3
fi

hisList=${ALIEN_JDL_HISTOLIST-$hisList}

#
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/AliTPCcalibAlignInterpolationMacro.C"

#
# we need macro locally to compile it
locMacro=$(basename "$macroName")
if [[ ! -f "$locMacro" ]] ; then cp $macroName ./ ; fi
#
#
for ihis in ${hisList[@]}
do 
    alilog_info "BEGIN: Processing map $ihis (${treeNames[$ihis]}) for time bin $mapStartTime : $mapStopTime in run $runNumber"
    #
    alilog_info "BEGIN: Projecting residual tree"
    time aliroot -b -q  $inclMacro ${locMacro}+\(1,$ihis,0\) >& projection_His${ihis}_Time${mapStartTime}.log
    alilog_info "END: Projecting residual tree"
    if [[ -f syswatch.log ]] ; then mv syswatch.log syswatch_His${ihis}_Time${mapStartTime}.log ; fi
    #
    alilog_info "BEGIN: Creating distortion map"
    time aliroot -b -q  $inclMacro ${locMacro}+\(2,$ihis,0\) >& map_His${ihis}_Time${mapStartTime}.log
    alilog_info "END: Creating distortion map"
    if [[ -f syswatch.log ]] ; then mv syswatch.log syswatch_Map${ihis}_Time${mapStartTime}.log ; fi
    #
    alilog_info "BEGIN: Creating NDLocal parameterizations"
    export inputFile="ResidualMapFull_${ihis}.root"
    export inputTree=${treeNames[$ihis]}

    for side in  0 1; do
	export varTheta0=$[side-1]
	export varTheta1=$[side];
        #
	for (( sec=0; sec<18; sec+=$secStep )); do
	    export varSec0=$sec
	    export varSec1=$[sec+secStep]
	    alilog_info "BEGIN: Processing side ${side}, sectors ${varSec0} ${varSec1}"
	    time aliroot -b -q  $inclMacro ${locMacro}+\(3\) >& NDlocal__His${ihis}_Time${mapStartTime}_Side${side}_Sectors${varSec0}_${varSec1}.log
	    alilog_info "END: Processing side ${side}, sectors ${varSec0} ${varSec1}"
	    #
	done
    done
    alilog_info "END: Creating NDLocal parameterizations"
    alilog_info "END: Processing map $ihis (${treeNames[$ihis]}) for time bin $mapStartTime : $mapStopTime in run $runNumber"
    #
done


