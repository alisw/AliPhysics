#!/bin/bash

FILESTOCHECK="rec.log AliESDs.root Run*.ESD.tag.root"
REMOVEROOTONERROR=""

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

if [ -f OCDB.generating.job ]; then
    echo "* This was a special OCDB.root job for which I'll skip the rest of the validation" >> stdout

    mv stdout stdout.ocdb.log 2>/dev/null
    mv stderr stderr.ocdb.log 2>/dev/null

    if [ -f OCDB.root ]; then
        echo "* ODCB.root found" >> stdout.ocdb.log
        exit 0
    else
        echo "* Error: OCDB.root NOT found! Failing validation" >> stdout.ocdb.log
        echo "Validation error cause: OCDB.root not found" >> .alienValidation.trace
        exit 1
    fi
fi

error=0

if [ -f validation_error.message ]; then
    echo "*! Found an error message from the job: `cat validation_error.message`" >> stdout
    echo "Validation error cause: job exited with message: `cat validation_error.message`" >> .alienValidation.trace
    error=4
fi

fullLogList=""

for logfile in *.log stdout stderr; do
    if [ -f "$logfile" -a ${logfile:0:8} != "syswatch" ]; then
        fullLogList="$fullLogList $logfile"
    fi
done

echo "* Checking for errors in $fullLogList"

if grep -q -s "W-AliReconstruction::Run: No events passed trigger selection" $fullLogList; then
    echo "* Job is forcefully validated since no events passed trigger selection" >> stdout
else
    if [ -f "qa.log" ]; then
        FILESTOCHECK="$FILESTOCHECK QAresults.root EventStat_temp.root"
    fi
    
    MISSINGFILES=""
    
    for file in $FILESTOCHECK `cat validation_extrafiles.list 2>/dev/null`; do
        if [ ! -f "$file" ]; then
            if [ ! -z "$MISSINGFILES" ]; then
                MISSINGFILES="$MISSINGFILES,"
            fi
            
            MISSINGFILES="${MISSINGFILES}${file}"
        fi
    done
        
    if [ ! -z "$MISSINGFILES" ]; then
        error=1
        echo "Validation error cause: Required file(s) not found in the output: $MISSINGFILES" >> .alienValidation.trace
        echo "* Error: Required file(s) not found in the output: $MISSINGFILES" >> stdout
    fi
    
    cat <<EOF >grep.cnf
std::bad_alloc
Segmentation violation
Segmentation fault
Bus error
floating point exception
Killed
busy flag cleared
Cannot Build the PAR Archive
*** glibc detected ***
E-AliCDBGrid::PutEntry:
F-AliCDBGrid::
E-TAlienFile::ReadBuffer: The remote (removed) file is not open
EOF

    found=`grep -f grep.cnf -i -n $fullLogList | grep -v -m 1 "Particle is stuck"`
    
    if [ ! -z "$found" ]; then
        (
            echo "* Found an error message in the logs:"
            echo "$found"
            echo ""
        ) >> stdout
    
        (
            echo "Validation error cause: Standard error message(s) found in the logs:"
            echo "$found"
        ) >> .alienValidation.trace
        
        error=2
    fi

    cat <<EOF >grep.cnf
Abort
Break
EOF
    
    if [ $error -eq 0 ]; then
        found=`grep -f grep.cnf -w -n -m 1 $fullLogList`
    
        if [ ! -z "$found" ]; then
    	    (
                echo "* Found an error message in the logs:"
                echo "$found"
                echo ""
    	    ) >> stdout
            
            (
                echo "Validation error cause: Standard error keyword(s) found in the logs"
                echo "$found"
            ) >> .alienValidation.trace
    
    	    error=2
        fi
    fi
    
    rm grep.cnf
    
    if [ -f check.log ]; then
        if ! grep -q "check of ESD was successfull" check.log; then
            echo "*! The ESD was not successfully checked according to check.log" >> stdout
            
            echo "Validation error cause: The ESD was not successfully checked according to check.log" >> .alienValidation.trace
            
            error=3
        fi
    fi
fi

(
if [ $error -eq 0 ]; then
    echo "* ################   Job validated ####################"
else
    echo "*! ################   Job NOT validated, error code $error ################"
    
    if [ ! -z "$REMOVEROOTONERROR" ]; then
        echo "* ########## Removing all ROOT files from the local directory, leaving only the logs ###"
        rm -rf *.root
    fi
    
    echo "* ulimit settings"
    ulimit -a
    echo "* memory limits"
    free -m
    echo "* last dmesg lines"
    dmesg | tail -n 20
    echo "* the rest of the environment settings"
    env
    echo "* End of machine status dump"
fi
) >> stdout

mv stdout stdout.log
mv stderr stderr.log

#(
#    echo "* Full env is:"
#    env | sort
#) >> stdout.log

exit $error
