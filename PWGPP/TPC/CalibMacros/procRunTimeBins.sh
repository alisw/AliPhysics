#!/bin/bash
# 
# This is the scropt to process multiple time bins of already preprocessed 
# residuals calibration objec
#
# arguments: 
# 1) timebins: log file with 1 line per timebin: tmin tmax run
Usage() {
    echo "Usage: ${0##*/} <timebins.log> "
    exit
}

[[ $# -lt 1 ]] &&  Usage && exit

timebins=$1

source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh  
#
loadLibMacro="$ALICE_PHYSICS/PWGPP/CalibMacros/CPass1/LoadLibraries.C"
inclMacro="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/includeMacro.C"
macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/procResidData.C"
prepobj="alitpcdcalibres.root"

# check if the preprocessed object exists
if [[ ! -e $prepobj ]] ; then
    alilog_info "Error: preprocessed object $prepobj is not in working directory"
    exit
fi

# we need macro locally to compile it
locMacro=$(basename "$macroName")
[[ ! -f "$locMacro" ]] && cp $macroName ./ 
[[ ! -f "$locMacro" ]] && echo "did not find $locMacro" && exit

bins=`wc -l $timebins | cut -f1 -d' '`
alilog_info "found $bins timebin in $timebins"
count=0
curdir=`pwd`

while read -r line ; 
do
    count=$[count + 1]
    alilog_info "Processing line $count of $bins: $line" 
    vals=($line)
    if [[ ${#vals[*]} -eq 3 ]] ; then
	tmin=${vals[0]}
	tmax=${vals[1]}
	runNumber=${vals[2]}
	run=$(echo "$runNumber" | sed 's/^0*//')
	binDir="$curdir"/"$tmin"_"$tmax"_"$runNumber"
	[ ! -d $binDir  ] && mkdir $binDir
	alilog_info "Processing bin $tmin - $tmax of run $run to directory $binDir"
	cp ${prepobj} $binDir/
	cd ${binDir}	
	mode=2
	time aliroot -b -q  $inclMacro $loadLibMacro ${curdir}/${locMacro}+g\($mode,$run,$tmin,$tmax,\"\"\) >& out_${tmin}_${tmax}.log
	rm tmp*.root
	alilog_info "Done"
    else
	alilog_info "Expect 3 arguments per line got ${#vals[*]} from $line"	
    fi
done < $timebins

exit 0
