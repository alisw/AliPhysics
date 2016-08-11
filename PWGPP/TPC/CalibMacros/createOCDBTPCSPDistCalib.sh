#!/bin/bash
#
# Shell script to generate the OCDB TPC calibration for SP distortions from the
# Chebychev polynomials
#
# arguments: file (can be an xml collection, or a txt file), start run for
# OCDB entry, end run for OCDB entry (optional: if not specified, the same
# as the start run will be used --> single run validity!), OCDB folder where
# to upload the file (optional, the default is the current local directory)
# Extra argument is what to process:
# corr=xxx : create correction object (xxx is any nonempty word not equal to "0")
# dist=xxx : create distortion object (xxx is any nonempty word not equal to "0")
# if both "corr" and "dist" are absent then only "correction" is processed
#
# e.g.: ./CreateOCDBTPCSPDistCalib.sh inputFileList=fileList.xml startRun=245731 endRun=245780 ocdbStorage="raw://" [corr=1] [dist=1]
#

###########################################################################

Usage() {

    echo "Usage ${0##*/}: wrong argument list, requires minimum the definition of inputFileList and startRun, e.g.: inputFileList=fileList.xml startRun=245731 "
    return 0
}

parseConfig()
{
    echo ""
    echo "Parsing the arguments"
    args=("$@")
    for opt in "${args[@]}"; do
	[[ -z ${opt} ]] && continue
	[[ -n ${encodedSpaces} ]] && opt="$(decSpaces ${opt})"
	[[ "${opt}" =~ ^[[:space:]]*$ ]] && continue
	if [[ ! "${opt}" =~ .*=.* ]]; then
	    echo "badly formatted option \"${opt}\" should be: option=value, stopping..."
	    return 1
	fi
	local var="${opt%%=*}"
	local value="${opt#*=}"
	export ${var}="${value}"
	echo "${var}=${value}"
    done
    
    return 0
}

arguments=("$@")
echo "The arguments are: $@"
#echo "${arguments[*]}"
endRun=-1
targetOCDBDir="local://"`pwd`
parseConfig "$@"
echo ""
if [[ -z $inputFileList || -z $startRun ]]; then 
    Usage
    return 0
fi

# allowing JDL to overwrite the default folder where to store the calibration 
targetOCDBDir=${ALIEN_JDL_TARGETOCDBDIR-$targetOCDBDir}

# check if correction/distortion objects were explicitly requested
[[ -z $corr || "$corr" == "0" ]] && corr="0" || corr="1"
[[ -z $dist || "$dist" == "0" ]] && dist="0" || dist="1"

# if neither of corr or dist is defined, impose correction only
[[ "$corr" == "0" && "$dist" == "0" ]] && corr="1"


echo "inputFileList = $inputFileList"
echo "startRun      = $startRun"
echo "endRun        = $endRun"
echo "ocdbStorage   = $ocdbStorage"
echo "corr          = $corr"
echo "dist          = $dist"


startRun=$(echo "$startRun" | sed 's/^0*//')
endRun=$(echo "$endRun" | sed 's/^0*//')

macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/ProcessOutputCheb.C"
locMacro=$(basename "$macroName")
if [[ ! -f "$locMacro" ]] ; then cp $macroName ./ ; fi

[ -e ocdb.log ] && rm ocdb.log

if [ "$corr" == "1" ]; then
    time aliroot -b -q "${locMacro}(\"$inputFileList\", $startRun, $endRun, \"$targetOCDBDir\",1)" >> ocdb.log
fi
#
if [ "$dist" == "1" ]; then
    time aliroot -b -q "${locMacro}(\"$inputFileList\", $startRun, $endRun, \"$targetOCDBDir\",0)" >> ocdb.log
fi


												      
