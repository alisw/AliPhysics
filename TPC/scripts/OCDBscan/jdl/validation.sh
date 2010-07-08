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

echo "... testing if the suppossed output was created ">> stdout
if [ -f calib.log ] && [ -f calib*.root ] 
    then 
    echo "ok">> stdout
    echo "... testing for Segmentation violation">> stdout
    sv=`grep -i  "Segmentation violation" *.log`
    if [ "$sv" = "" ]
	then
	echo "ok">> stdout
	echo "... testing for Segmentation fault">> stdout
	sf=`grep -i  "Segmentation fault" *.log`
	if [ "$sf" = "" ]
	    then
	    echo "ok">> stdout
	    echo "... testing for Bus error">> stdout
	    be=`grep -i  "Bus error" *.log`
	    if [ "$be" = "" ]
		then
		echo "ok">> stdout
		echo "... testing for Break">> stdout
		ab=`grep -i "Break" *.log`
		if [ "$ab" = "" ]
		    then
		    echo "ok">> stdout
		    echo "... testing for Floating point exception">> stdout
		    fp=`grep -i  "Floating point exception" *.log`
		    if [ "$fp" = "" ]
			then
			echo "ok">> stdout
			echo "... testing for Killed">> stdout
			kl=`grep -i  "Killed" *.log`
			if [ "$kl" = "" ]
			    then
			    echo "ok">> stdout
			    echo "... testing for busy flag cleared">> stdout
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
    echo "* ########## Job not validated - no calib.log or calibSummary.root ###" >> stdout;
    echo "* ########## Removing all ROOT files from the local directory, leaving only the logs ###" >> stdout;
    rm -rf *.root
fi


if [ "$error" = "1" ] 
    then
    echo "* ################   Job not validated ################" >> stdout;
fi


echo "* ----------------------------------------------------*" >> stdout;
echo "*******************************************************" >> stdout;


sleep 15;
cd -

cp stdout stdout.log
if [ -f stderr ]
    then
    cp stderr stderr.log
fi

exit $error
