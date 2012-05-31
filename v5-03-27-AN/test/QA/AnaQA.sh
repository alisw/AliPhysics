#!/bin/sh

# AnaQA.sh
# 
#
# Created by schutz on 30/09/08.
# Copyright 2008 CERN. All rights reserved.
macroname="AnaQA"
validateout=`dirname $0`
validatetime=`date`
validated="0";
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout;
echo "* AliRoot QA Validation Script V1.0                   *" >> stdout;
echo "* Time:    $validatetime " >> stdout;
echo "* Dir:     $validateout" >> stdout;
echo "* Workdir: $validateworkdir" >> stdout;
echo "* ----------------------------------------------------*" >> stdout;
detectorlist="ITS TPC TRD TOF PHOS HMPID EMCAL FMD ZDC T0 VZERO PMD ACORDE Global"
if [ $# -eq 0 ] ; then
 echo "syntax: AnaQA.sh #runNumber"
 exit 1
fi
if [ ! -f  ${macroname}.C ] ; then 
 echo "* ########## Job not validated - no validation macro (${macroname}.C)  ###" >> stdout;
 exit 2
fi 
run=$1
logfile=${macroname}_${run}.log
if [ -e $logfile ] ; then
 rm $logfile
fi 
aliroot -b > $logfile <<EOF
.L ${macroname}.C+
${macroname}($run)
EOF
if [ ! -f  $logfile ] ; then 
 echo "* ########## Job not validated - no log file ($logfile)  ###" >> stdout;
 exit 2
fi 
let error=0
declare -a array
for pb in `grep -i "Problem signalled" $logfile | awk '{print $1"_"$2"_"$NF}'` ; do
 for det in detectorlist ; do
  array[$error]=`echo $pb | awk '{ split($0, a, "_"); print a[1]" "a[3]" in "a[2]}'`
  ((error++))
 done 
done  
if [ "$error" = "0" ] ; then 
 echo "* ----------------   Job Validated  ------------------*" >> stdout;
else 
 errors=${#array[@]}
 echo "* ########## Job not validated - number of errors: $errors ###" >> stdout;
 for (( i=0;i<$errors;i++ )); do
  echo $i-- ${array[${i}]} >> stdout;
 done
fi
echo "* ----------------------------------------------------*" >> stdout;
echo "*******************************************************" >> stdout;
