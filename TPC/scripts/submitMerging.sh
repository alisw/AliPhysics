# submit merging of the calibration train
# argument 1 - runlist
# argument 2 - batch queues

runlist=$1
bqueues=$2
mydir=`pwd`
for adir in `cat $runlist`; do
    myvar=0;
    cd $mydir
    cd $adir
    echo Run $adir 
    nesd=`wc -l < esd.txt.Good`
    if [ $nesd -gt 0 ] ; then
	rm -f CalibObjects*.root
	find  `pwd`/*_*  | grep CalibObjectsTrain1.root > mergelistTrain1.txt
	find  `pwd`/*_*  | grep CalibObjectsTrain2.root > mergelistTrain2.txt
	nfiles=`cat mergelistTrain1.txt  | grep -c .root`
	if [ $nfiles -gt 0 ] ; then
	    bsub -q $bqueues -oo outm1_$myvar.log aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergelistTrain1.txt\"\)
	    bsub -q $bqueues  -oo outm2_$myvar.log aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergelistTrain2.txt\"\)
	    echo Run $adir  Nfiles=$nfiles
        else
            echo Run $adir Calib Missing
        fi;
    fi; 
    cd $mydir
done
 
