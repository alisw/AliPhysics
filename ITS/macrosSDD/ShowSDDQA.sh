#!/bin/bash
RUN='123456'
PERIOD='LHC10a'
PASS='pass1'
YEAR='2010'
ALICE_ITS='$ALICE_ROOT/ITS'
echo "Run Number   :[${RUN}]"
read
if [ "$REPLY" != "" ]; then
export RUN=$REPLY
echo "Run   $RUN"
fi
echo "Period        :[${PERIOD}]"
read
if [ "$REPLY" != "" ]; then 
export PERIOD=$REPLY
echo "Period  $PERIOD"
fi
echo "Pass           : [${PASS}]"
read
if [ "$REPLY" != "" ]; then
export PASS=$REPLY
echo "Pass       $PASS "
fi
echo "Year      :[${YEAR}]"
read
if [ "$REPLY" != "" ]; then
export YEAR=$REPLY
echo "Year    $YEAR"
fi
if [ ls -l "run$RUN" >/dev/null 2>&1 ]; then
echo "directory run$RUN exists "
else
mkdir "run$RUN"
fi
cd "run$RUN"
if [ ls -l $PASS > /dev/null 2>&1 ]; then
echo "directory $PASS exixsts"
else
mkdir $PASS
cd $PASS
fi
time aliroot >>merge.log 2>&1 <<EOI
.x $ALICE_ITS/ReadQASDD.C($RUN,$YEAR,"${PERIOD}","${PASS}" ); 
.q
EOI
time aliroot  >> plot.log 2>&1 <<EOI
.x $ALICE_ITS/PlotQASDD.C("File.QA.${YEAR}.${PERIOD}.${PASS}.Run.${RUN}.root");
EOI

if [ls -l "images" >/dev/null 2>&1 ]; then
echo "directory images exists"
else
mkdir images
fi
mv *.ps *.eps images/.
cd images
for i in  *.ps;
do
gv $i & 
sleep 2
done
echo "Plots Done!!"

