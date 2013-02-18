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
echo "* Library path: $LD_LIBRARY_PATH" >> stdout;
echo "* ----------------------------------------------------*" >> stdout;
echo "* Path: $PATH" >> stdout;
echo "* ----------------------------------------------------*" >> stdout;
ls -la ./ >> stdout;
echo "* ----------------------------------------------------*" >> stdout;

##################################################
#if [ -f rec.log ] && [ -f sim.log ] && [ -f check.log ] && [ -f tag.log ] && [ -f aod.log ] && [ -f *ESD.tag.root ]
#if [ -f rec.log ] && [ -f sim.log ] && [ -f check.log ] && [ -f tag.log ] && [ -f aod.log ] && [ -f *ESD.tag.root ] && [ -f AnalysisResults.root ]
if [ -f rec.log ] && [ -f sim.log ] && [ -f checkesd.log ] && [ -f checkaod.log ] && [ -f aod.log ]
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
            ab=`grep -i -w "Abort" *.log`
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
                            ch=`grep -i "check of ESD was successfull" checkesd.log`
                            if [ "$ch" = "" ]
                            then
                                echo "* #  The ESD was not successfully checked   *"  >>stdout;
                            else
                                ao=`grep -i "check of AOD was successfull" checkaod.log`
                                if [ "$ao" = "" ]
                                then
                                    echo "* #  The AOD was not successfully checked   *"  >>stdout;
                                else
                                    echo "* ----------------   Job Validated  ------------------*" >> stdout;
                                    error="0";
				fi
                            fi
                        fi
                    fi
                fi
            fi
        fi
    fi
fi
else
    echo "* ########## Job not validated - no rec.log or sim.log or checkaod.log or checkesd.log ###" >> stdout;
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
exit $error

