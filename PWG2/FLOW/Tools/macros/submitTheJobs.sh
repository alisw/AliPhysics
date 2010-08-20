#!/bin/sh

#### this is user settable #########
#how many files per job
n=4

#file with the list for each job
export datalistfile="files.txt"

#where do we get our data from?
find /data/alice3/jthomas/testData/ -name "*.root" > $datalistfile

########## DO NOT TOUCH THIS SECTION ###########################################
################################################################################
#script to be run
export scripttorun="theNodeScript.sh"

queue="stbcq"
[[ $# -eq 1 ]] && queue=$1
export listoffiles=""
i=0
for x in `cat $datalistfile`
do
  listoffiles="$listoffiles $x"
  if [ $(($i%$n)) -eq 0 ]
  then
    echo $listoffiles
    echo ""
    qsub -q ${queue} -v "listoffiles,datalistfile" $scripttorun
    listoffiles=""
  fi
  ((i++))
done
################################################################################
