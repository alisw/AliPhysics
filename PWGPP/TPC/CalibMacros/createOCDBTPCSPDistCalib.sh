#!/bin/bash
#
# Shell script to generate the OCDB TPC calibration for SP distortions from the
# Chebychev polynomials
#
# arguments: file (can be an xml collection, or a txt file), start run for
# OCDB entry, end run for OCDB entry (optional: if not specified, the same
# as the start run will be used --> single run validity!), OCDB folder where
# to upload the file (optional, the default is the current local directory)
#
# e.g.: ./CreateOCDBTPCSPDistCalib.sh inputFileList=fileList.xml startRun=245731 endRun=245780 ocdbStorage="raw://"
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

echo "inputFileList = $inputFileList"
echo "startRun      = $startRun"
echo "endRun        = $endRun"
echo "ocdbStorage   = $ocdbStorage"

startRun=$(echo "$startRun" | sed 's/^0*//')
endRun=$(echo "$endRun" | sed 's/^0*//')

macroName="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/ProcessOutputCheb.C"
locMacro=$(basename "$macroName")
if [[ ! -f "$locMacro" ]] ; then cp $macroName ./ ; fi

time aliroot -b -q "${locMacro}(\"$inputFileList\", $startRun, $endRun, \"$targetOCDBDir\")" >& ocdb.log


												      
