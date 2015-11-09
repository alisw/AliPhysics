#!/bin/bash

FILESTOCHECK="rec_Barrel.log rec_Outer.log calib.log AliESDs_Barrel.root CalibObjects.root"

validateout=`dirname $0`

if [ -z "$validateout" ]; then
    validateout="."
fi

cd "$validateout"
validateworkdir=`pwd`

(
echo "* *****************************************************"
echo "* AliRoot Validation Script V2.0                      *"
echo "* Time:    `date`"
echo "* Dir:     $validateout"
echo "* Workdir: $validateworkdir"
echo "* ----------------------------------------------------*"

for subdir in . Barrel OuterDet; do
    if [ -d "$subdir" ]; then
        echo "Listing $subdir" 
        ls -lA "$subdir/"
        echo ""
    fi
done

echo "* ----------------------------------------------------*"
) >> stdout

if [ -f OCDB.generating.job ]; then
    echo "* This was a special OCDB.root job for which I'll skip the rest of the validation" >> stdout

    mv stdout stdout.ocdb.log 2>/dev/null
    mv stderr stderr.ocdb.log 2>/dev/null
    mv rec.log stdrec.ocdb.log 2>/dev/null

    if [ -f OCDB.root ]; then
        echo "* ODCB.root found" >> stdout.ocdb.log
        exit 0
    else
        echo "* Error: OCDB.root NOT found! Failing validation" >> stdout.ocdb.log
        exit 1
    fi
fi

error=0

for file in $FILESTOCHECK; do
    if [ ! -f "$file" ]; then
        error=1
        echo "* Error: Required file $file not found in the output" >> stdout
    fi
done

for message in "std::bad_alloc" "Segmentation violation" "Segmentation fault" "Bus error" "Break" "Floating point exception" "Killed" "busy flag cleared"; do
    found=`grep -i -n -m 1 "$message" *.log`
    if [ ! -z "$found" ]; then
        (
        echo "* Found error message '$message' in the logs:"
        echo "$found"
        echo ""
        ) >> stdout

        error=2
    fi
done

(
if [ $error -eq 0 ]; then
    echo "* ################   Job validated ####################"
else
    echo "* ################   Job NOT validated, error code $error ################"
#    echo "* ########## Removing all ROOT files from the local directory, leaving only the logs ###"
#    rm -rf *.root
fi
) >> stdout

mv stdout stdout.log
mv stderr stderr.log

exit $error
