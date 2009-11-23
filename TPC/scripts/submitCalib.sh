#
# Submit calibration in predefined workspace
# argument 1 - runlist 
# argument 2 - batch queues
# arument  3 - step
# Assumptions:
# 1. Workspace defined before using the script $ALICE_ROOT/TPC/scripts/makeWorkspace.sh
# 2. The data were filtered before using the $ALICE_ROOT/TPC/scripts/filterCorrupted.sh

runlist=$1
bqueues=$2
step=$3
echo runlist"      "$runlist
echo queues"       "$bqueues
echo step"         "$step
echo
mydir=`pwd`
for adir in `cat $runlist`; do
    cd $mydir
    #
    #remove old data
    rm -f *.log
    rm -f core*
    rm -rf $mydir/$adir/core*
    rm -rf $mydir/$adir/*/core*
    rm -rf $mydir/$adir/*.root
    rm -rf $mydir/$adir/*.log
    rm -rf $mydir/$adir/*_*
    #end of remove old data
    #
    up=`cat  $adir/esd.txt.Good | grep -c .root`
    if [ $up -gt 0 ] ; then
	myvar=0;
	cd $adir
	echo SUBMITING DIRECTORY $adir
	cp ../ConfigOCDB.C .
	cp ../CalibrateTPC.C .
	while [ $myvar -le ${up} ] ; do
        bsub -q $bqueues  -oo out$myvar.log -eo err$myvar.log $ALICE_ROOT/TPC/scripts/submitCalibJob.sh $myvar $(($myvar+step))  `pwd`/esd.txt.Good $adir; 
	myvar=$(( $myvar + $step )) ; 
	echo $myvar ; 
    done;
fi;
done;
