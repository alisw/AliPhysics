#!/bin/bash
currentDir=`pwd`
PACKAGES="ESD ANALYSIS AnalysisCheck"
DETECTORS="PHOS EMCal PMD HMPID T0 MUON TOF VZERO"
LOGFILE=`echo $0 | sed -e 's/sh/log/'`
function testpack()
{
 test=`cat CheckESD.log | grep "$2 done"` 
 if [ "$test" = "" ]; then 
    echo --ERROR $1 $2
    rv=1 
 else 
    echo ++OK $1 $2
    rv=0 
 fi
 return $rv
}
function testdet()
{
 test=`cat CheckESD.log | grep "$2 Summary Report: OK"`
 if [ "$test" = "" ]; then
    echo --ERROR $1 $2 
    rv=1
 else 
    echo ++OK $1 $2
    rv=0
 fi
 return $rv
}
#Start
if [ -e $LOGFILE ]; then
 rm $LOGFILE
fi
echo $0 LOG > $LOGFILE
# make the par files
cd $ALICE_ROOT
for pack in $PACKAGES; do
 make $pack.par            >> $currentDir/$LOGFILE
done
# copy the par file to the working directory
cd $currentDir
for pack in $PACKAGES; do
 rm -fr $pack* 
 mv $ALICE_ROOT/$pack.par .
done
cp $ALICE_ROOT/ESDCheck/ana.C .
# run root
if [ ! -e "AliESDs.root" ]; then 
 echo File AliESDs.root not found >> $LOGFILE
 exit 1 
fi
root -b -q ana.C >> $LOGFILE 2>&1
#test function for parsing log file
# parse the log file
error=0
for pack in $PACKAGES; do
 testpack creating $pack.par  
 testpack loading lib$pack 
 let "error +=$?"
done

for det in $DETECTORS; do
 testdet Checking $det
 let "error +=$?"
done 
if [ "$error" > "0" ]; then 
 echo --------------- $error errors signaled 
else 
 echo +++++++++++++++ All OK
fi 
exit $error
 
