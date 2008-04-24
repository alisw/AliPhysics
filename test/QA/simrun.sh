#!/bin/csh
if ($#argv < 1) then
 echo "usage simrun.sh RunNumber"
 exit()
endif
cd $WORK/QATest
rm -f *.root *.C *.log
ln -si $ALICE_ROOT/test/QA/Config.C Config.C
ln -si $ALICE_ROOT/test/QA/sim.C sim.C
ln -si $ALICE_ROOT/test/QA/simqa.C simqa.C
ln -si $ALICE_ROOT/test/QA/rec.C rec.C
ln -si $ALICE_ROOT/test/QA/recqa.C recqa.C
ln -si $ALICE_ROOT/GRP/Geometry/Data/Run0_999999999_v0_s0.root geometry.root
root -b -q $ALICE_ROOT/test/QA/simrun.C --run $1
cd $WORK/QATest/data
rm -f *.root
cp  $ALICE_ROOT/test/QA/rawqa.C .
ln -si $ALICE_ROOT/GRP/Geometry/Data/Run0_999999999_v0_s0.root geometry.root
alienaliroot -b > rawqa.log 2>&1 << EOF
.x  $ALICE_ROOT/test/QA/rootlogon.C
.L rawqa.C++
rawqa(21950, 10)
EOF
rm -f rawqa.C
exit
