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
    up=`cat  $adir/esd.txt.Good | grep -c .root`
    if [ $up -gt 1 ] ; then
	myvar=0;
	cd $adir
	echo SUBMITING DIRECTORY $adir
	cp ../ConfigOCDB.C .
	cp ../CalibrateTPC.C .
	rm -rf *_*
	rm *.log
	rm -rf */V3/
	while [ $myvar -le ${up} ] ; do
        bsub -q $bqueues  -oo out$myvar.log -eo err$myvar.log $ALICE_ROOT/TPC/scripts/submitCalibJob.sh $myvar $(($myvar+step))  `pwd`/esd.txt.Good $adir; 
	myvar=$(( $myvar + $step )) ; 
	echo $myvar ; 
    done;
fi;
done;
