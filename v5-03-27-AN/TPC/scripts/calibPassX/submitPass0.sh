#!/bin/bash
#
# Submit jobs for reconstruction/calibration
# Test - Equivalent of pass0/passX on the grid
#

# Parameters:
# 1 - queue name
# 2 - run number
# 3 - number of events to be reconstructed
# 4 - run.list
# Example: submitPass0.sh run.list "alice-t3_8h" 1000 

runList=$1
qName=$2
nEvents=$3
workDir=`pwd`
for arun in `cat $runList`; do
    echo 
    echo Run number $arun
    echo 
    #
    cd $workDir/$arun
    cp $workDir/*.C .
    cp $workDir/*.sh .
    cat ../../lists/raw.list | grep $arun >raw.list    
    #
    #
    for afile in `cat raw.list | grep -v tag`; do 
	echo $afile; 
	cdir0=`basename $afile | sed s_.root__`
	mkdir $cdir0
	cd $cdir0
	cp $workDir/$arun/*.C .
	cp $workDir/$arun/*.sh .
        echo
        echo
        pwd
	echo bsub -q $qName -o out${cdir0}.log   runPassXJob.sh $afile  $nEvents  $arun
	bsub      -q $qName -o out${cdir0}.log   runPassXJob.sh $afile  $nEvents  $arun 
	cd ../
        cd $workDir/$arun
    done;
    cd $workDir
done



