#!/bin/bash

for rcfile in validation.rc validation_merge.rc; do
    if [ -f $rcfile ]; then
        source $rcfile
    fi
done

validateout=`dirname $0`

if [ -z "$validateout" ]; then
    validateout="."
fi

cd "$validateout"
validateworkdir=`pwd`

(
    echo "* *****************************************************"
    echo "* AliRoot Validation Script V2.3                      *"
    echo "* Time:    `date`"
    echo "* Dir:     $validateout"
    echo "* Workdir: $validateworkdir"
    echo "* PATH: $PATH"
    echo "* LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
    echo "* ----------------------------------------------------*"
    
    for subdir in . Barrel OuterDet; do
	if [ -d "$subdir" ]; then
            echo "Listing $subdir" 
            ls -lA "$subdir/"
            echo ""
	fi
    done
    
    echo "df -h ."
    df -h .
    
    echo "* ----------------------------------------------------*"
) >> stdout

error=0

if [ -f check.log ]; then
    if ! grep -q "check of ESD was successfull" check.log; then
	echo "*! The ESD was not successfully checked according to check.log" >> stdout
	
	echo "Validation error cause: The ESD was not successfully checked according to check.log" >> .alienValidation.trace
	
	error=4
    fi
fi

if ! [ -f jstiller_AliESDs.root ] ; then
    error=1
    echo "Output file jstiller_AliESDs.root not found. Job FAILED !" >> stdout
    echo "Output file jstiller_AliESDs.root not found. Job FAILED !" >> stderr
fi
if ! [ -f AliAOD.root ] ; then
    error=2
    echo "Output file AliAOD.root not found. Job FAILED !" >> stdout
    echo "Output file AliAOD.root not found. Job FAILED !" >> stderr
fi
if ! [ -f AliAOD.VertexingHF.root ] ; then
    error=3
    echo "Output file AliAOD.VertexingHF.root not found. Job FAILED !" >> stdout
    echo "Output file AliAOD.VertexingHF.root not found. Job FAILED !" >> stderr
fi

(
    if [ $error -eq 0 ]; then
	echo "* ################   Job validated ####################"
    else
	echo "*! ################   Job NOT validated, error code $error ################"
    fi
) >> stdout

mv stdout stdout.log
mv stderr stderr.log


exit $error
