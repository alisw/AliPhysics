#!/bin/sh -f 
if ! [ "$1" ]
then
 echo "Usage: simrun.sh RunNumber";
 exit;
fi

# Working directory
if [ ${WORK}"x" == "x" ]
then 
 export WORK=./;
fi
cd ${WORK};

# Version
if [ ${VERSION}"x" == "x" ]
then 
 export VERSION=UNKNOWN;
fi
echo $WORK $VERSION

# Test directory
if ! [ -e QATest/$VERSION ]
then 
 mkdir -p QATest/$VERSION;
fi
cd QATest/$VERSION;

# Cleanup and (re)link
rm -rf *.root *.C *.log data;
ln -sf $ALICE_ROOT/test/QA/Config.C Config.C;
ln -sf $ALICE_ROOT/test/QA/sim.C sim.C;
ln -sf $ALICE_ROOT/test/QA/simqa.C simqa.C;
ln -sf $ALICE_ROOT/test/QA/rec.C rec.C;
ln -sf $ALICE_ROOT/test/QA/recqa.C recqa.C;

# Process MC
root -b -q $ALICE_ROOT/test/QA/simrun.C --run $1;

# Directory for RAW data 
if [ ! -e data ]
then 
 mkdir data;
fi
cd data;

# Cleanup and (re)link
rm -rf *.root *.C *.log;
ln -sf ../raw.root;
ln -sf $ALICE_ROOT/test/QA/recraw.C recraw.C;
cp -f $ALICE_ROOT/test/QA/rawqa.C .;

# Process RAW
aliroot -b -q recraw.C  > recraw.log ;
if ! [ ${GSHELL_ROOT}"x" == "x" ] # Check that the AliEn API is available
then
aliroot -b > rawqa.log << EOF
.x  $ALICE_ROOT/test/QA/rootlogon.C
.L rawqa.C++
rawqa($1, 10)
EOF
fi
exit;
