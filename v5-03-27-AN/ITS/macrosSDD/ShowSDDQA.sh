#!/bin/bash
RUN='123456'
PERIOD='LHC10a'
PASS='pass1'
YEAR='2010'
ALICE_ITS='$ALICE_ROOT/ITS'
TMPPLACE='/tmp'
TMPFOLDER='1'
EXECFOLDER='$HOME/macroQAshifter'
MAXFILES='300'
echo "Run Number   :[${RUN}]"
read
if [ "$REPLY" != "" ]; then
RUN=$REPLY
echo "Run   $RUN"
fi
echo "Period        :[${PERIOD}]"
read
if [ "$REPLY" != "" ]; then 
PERIOD=$REPLY
echo "Period  $PERIOD"
fi
echo "Pass           : [${PASS}]"
read
if [ "$REPLY" != "" ]; then
PASS=$REPLY
echo "Pass       $PASS "
fi
echo "Year      :[${YEAR}]"
read
if [ "$REPLY" != "" ]; then
YEAR=$REPLY
echo "Year    $YEAR"
fi
echo "folder with macros     :[${EXECFOLDER}]"
read
if [ "$REPLY" != "" ]; then
EXECFOLDER=$REPLY
echo "Folder:    $EXECFOLDER"
fi
echo "local or lxplus (1=local 2=lxplus)   :[${TMPFOLDER}]"
read
if [ "$REPLY" != "" ]; then
TMPFOLDER=$REPLY
fi
if [ "$TMPFOLDER" == "1" ]; then
TMPPLACE='/tmp'
else
TMPPLACE='/tmp/$USERNAME'
fi
GOOD=130
echo "Max number of files (Insert a number >0)   :[${MAXFILES}]"
read
#case $REPLY in
#	*[a-zA-Z]*|*[!0-9]*|*[-]*)
#	                            echo "Wrong: it is not a number > 0";;   
#			  *[0-9]*)
#				   MAXFILES=$REPLY
#				   echo "Max number of files:    $MAXFILES" 
# 				   export GOOD=0
#			;;
#esac	
if [ "$REPLY" != "" ]; then
MAXFILES=$REPLY
echo "Max number of files:    $MAXFILES"
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
time aliroot -l <<EOI|tee merge$RUN.log
EOF
.x $EXECFOLDER/ReadQASDD.C($RUN,$YEAR,"${PERIOD}","${PASS}",$MAXFILES); 
.q
EOI
time aliroot -l <<EOI|tee plot$RUN.log
.x $EXECFOLDER/PlotQASDD.C("File.QA.${YEAR}.${PERIOD}.${PASS}.Run.${RUN}.root");
.q
EOI
rm File.QA.${YEAR}.${PERIOD}.${PASS}.Run.${RUN}.root
if [ls -l "images" >/dev/null 2>&1 ]; then
echo "directory images exists"
else
mkdir images
fi
mv *.ps images/.
cd images
for i in  *.ps;
do
gv $i & 
sleep 2
done
echo "Plots Done!!"
rm -rf $TMPPLACE/*.root
cd ../../../
