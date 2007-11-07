#!/bin/sh
##################################################
validateout=`dirname $0`
validatetime=`date`
validated="0";
error=1
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout;
echo "* AliRoot Validation Script V1.0                      *" >> stdout;
echo "* Time:    $validatetime " >> stdout;
echo "* Dir:     $validateout" >> stdout;
echo "* Workdir: $validateworkdir" >> stdout;
echo "* ----------------------------------------------------*" >> stdout;
ls -la ./ >> stdout;
echo "* ----------------------------------------------------*" >> stdout;

##################################################
if [ -f job.log ] 
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
	    ab=`grep -i "Abort" *.log`
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
			    echo "* ----------------   Job Not Validated  ------------------*" >> stdout;			    
			fi
		    fi;
		fi
	    fi
	fi
    fi
fi
else
    echo "* ########## Job not validated - no job.log  ###" >> stdout;
fi
if [ "$error" = "1" ] 
    then
    echo "* ################   Job not validated ################" >> stdout;
fi
echo "* ----------------------------------------------------*" >> stdout;
echo "*******************************************************" >> stdout;
sleep 15;
cd -
echo $error
exit $error
