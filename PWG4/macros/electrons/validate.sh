#!/bin/bash
##################################################
validateout=`dirname $0`
validatetime=`date`
validated="0";
error=0
if [ -z $validateout ]
then
    validateout="."
fi

cd $validateout;
validateworkdir=`pwd`;

echo "*******************************************************" >> stdout
echo "* Automatically generated validation script           *" >> stdout

echo "* Time:    $validatetime " >> stdout
echo "* Dir:     $validateout" >> stdout
echo "* Workdir: $validateworkdir" >> stdout
echo "* ----------------------------------------------------*" >> stdout
ls -la ./ >> stdout
echo "* ----------------------------------------------------*" >> stdout

##################################################
if ! [ -f histos.root ] ; then
   error=1
   echo "Output file(s) not found. Job FAILED !" >> stdout
   echo "Output file(s) not found. Job FAILED !" >> stderr
fi
if [ $error = 0 ] ; then
   echo "* ----------------   Job Validated  ------------------*" >> stdout
fi
echo "* ----------------------------------------------------*" >> stdout
echo "*******************************************************" >> stdout
cd -
exit $error
