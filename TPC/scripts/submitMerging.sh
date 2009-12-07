# submit merging of the calibration train
# argument 1 - runlist
# argument 2 - batch queues
# argument 3 - group

runlist=$1
bqueues=$2 
bgroup=$3

if [ -z $bgroup ] ; then
  bgroup=/merge
  bgadd $bgroup 
fi;

mydir=`pwd`
echo bqueues $bqueues
for adir in `cat $runlist`; do
    myvar=0;
    cd $mydir
    cd $adir
    echo Run $adir 
    nesd=`wc -l < esd.txt.Good`
    if [ $nesd -gt 0 ] ; then
	rm -f CalibObjects*.root
	find  `pwd`/*_*  | grep CalibObjectsTrain1.root | grep -v lxb  > mergelistTrain1.txt
	find  `pwd`/*_*  | grep CalibObjectsTrain2.root | grep -v lxb > mergelistTrain2.txt
	nfiles=`cat mergelistTrain1.txt  | grep -c .root`
	if [ $nfiles -gt 0 ] ; then
	    bsub -q $bqueues -g $bgroup -oo outm1_$myvar.log aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"mergelistTrain1.txt\"\)
	    bsub -q $bqueues -g $bgroup -oo outm2_$myvar.log aliroot $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"mergelistTrain2.txt\"\)
	    echo Run $adir  Nfiles=$nfiles
        else
            echo Run $adir Calib Missing
        fi;
    fi; 
    cd $mydir
done
 
