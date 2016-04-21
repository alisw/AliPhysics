#!/bin/bash

FILESTOCHECK="voxelResTree.root"
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

for subdir in . ; do
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

found=`grep -f grep.cnf -i -n $fullLogList`

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
