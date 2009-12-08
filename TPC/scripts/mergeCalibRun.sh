#
# merge calib data for given list
# input - run list
# Usage:
# $ALICE_ROOT/TPC/scripts/mergeCalibRun.sh run.list

runlist=$1
mydir=`pwd`
mkdir merge$runlist
cd  merge$runlist
rm train1.txt
rm train2.txt
touch train1.txt
touch train2.txt
for adir in `cat ../$runlist`; do
ls $mydir/$adir/CalibObjectsTrain1.root >> train1.txt;
ls $mydir/$adir/CalibObjectsTrain2.root >> train2.txt;
done;
aliroot -b -q  $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain1.root\",10000,\"train1.txt\"\)
aliroot -b -q  $ALICE_ROOT/TPC/macros/CalibFileMerger.C+\(\"CalibObjectsTrain2.root\",10000,\"train2.txt\"\)
cd ..
