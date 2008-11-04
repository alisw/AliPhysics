#!/bin/csh -f 
if ($#argv < 1) then
 echo "usage simrun.sh RunNumber"
 exit()
endif
echo $WORK $VERSION
if ( ! -e $WORK ) then 
 setenv WORK ./
endif
cd $WORK
if ( $VERSION == "" ) then 
 setenv VERSION "UNKNOWN"
endif
if ( ! -e QATest/$VERSION ) then 
 mkdir QATest/$VERSION
endif    
cd QATest/$VERSION
rm -Rf DB* *.root *.C *.log data/*
ln -si $ALICE_ROOT/test/QA/Config.C Config.C
ln -si $ALICE_ROOT/test/QA/sim.C sim.C
ln -si $ALICE_ROOT/test/QA/simqa.C simqa.C
ln -si $ALICE_ROOT/test/QA/rec.C rec.C
ln -si $ALICE_ROOT/test/QA/recqa.C recqa.C
ln -si $ALICE_ROOT/test/QA/DB.tgz DB.tgz
root -b -q $ALICE_ROOT/test/QA/simrun.C --run $1
if ( ! -e data ) then 
 mkdir data
endif
cd data
#ln -s ../geometry.root
ln -s ../raw.root
ln -s ../DB 
ln -si $ALICE_ROOT/test/QA/recraw.C recraw.C
aliroot -b -q recraw.C  > recraw.log 
cp  $ALICE_ROOT/test/QA/rawqa.C .
aliroot -b > rawqa.log << EOF
.x  $ALICE_ROOT/test/QA/rootlogon.C
.L rawqa.C++
rawqa($1, 10)
EOF
rm -f rawqa.C
exit
