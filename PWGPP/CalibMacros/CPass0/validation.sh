#!/bin/bash

validateout=`dirname $0`
validated="0";
error=1

if [ -z "$validateout" ]; then
    validateout="."
fi

cd "$validateout"
validateworkdir=`pwd`

(
echo "*******************************************************"
echo "* AliRoot Validation Script V2.0                      *"
echo "* Time:    `date`"
echo "* Dir:     $validateout"
echo "* Workdir: $validateworkdir"
echo "* ----------------------------------------------------*"
ls -lA .
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
        echo "* OCDB.root NOT found! Failing validation" >> stdout.ocdb.log
        exit 1
    fi
fi

cp stdout stdout.log
cp stderr stderr.log

##################################################
if [ -f rec.log ] && [ -f calib.log ] && [ -f AliESDs.root ] && [ -f CalibObjects.root ]  
then 
sv=`grep -i  "Segmentation violation" *.log`
if [ "$sv" = "" ]
    then
    sf=`grep -i  "Segmentation fault" *.log`
    if [ "$sf" = "" ]
        then
        be=`grep -i  "Bus error" *.log`
        if [ "$be" = "" ]
	    then
	    ab=`grep -i "Break" *.log`
	    if [ "$ab" = "" ]
	        then
	        fp=`grep -i  "Floating point exception" *.log`
	        if [ "$fp" = "" ]
		    then
		    kl=`grep -i  "Killed" *.log`
		    if [ "$kl" = "" ]
		        then
		        bf=`grep -i "busy flag cleared" *.log`
		        if [ "$bf" = "" ]
                            then
			       echo "* ----------------   Job Validated  ------------------*" >> stdout;
			       error="0";
                            else
                               echo "* #             Check Macro failed !                  #" >> stdout;
                        fi
		    fi
	        fi
            fi
        fi
    fi
fi
else
    echo "* ########## Job not validated - no rec.log or calib.log or AliESDs.root  && AliESDfriends.root ###" >> stdout;
    echo "* ########## Removing all ROOT files from the local directory, leaving only the logs ###" >> stdout;
    rm -rf *.root
fi
if [ "$error" = "1" ] 
    then
    echo "* ################   Job not validated ################" >> stdout;
fi
echo "* ----------------------------------------------------*" >> stdout;
echo "*******************************************************" >> stdout;

exit $error
