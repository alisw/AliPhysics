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
FILENAME='10000137137031.300.root'
FULLNAME='$PWD/$FILENAME'
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
echo "FileName (if it is a LOCAL file, please insert the full path)      :  [${FILENAME}]"
read
if [ "$REPLY" != "" ]; then
FILENAME=$REPLY
fi
echo "FileName    ${FILENAME}"
echo "folder with macros     :[${EXECFOLDER}]"
read
if [ "$REPLY" != "" ]; then
EXECFOLDER=$REPLY
echo "Folder:    $EXECFOLDER"
fi
echo "local file or alienfile (1=local 2=alien)   :[${TMPFOLDER}]"
read
if [ "$REPLY" != "" ]; then
TMPFOLDER=$REPLY
fi
if [ "$TMPFOLDER" == "1" ]; then
FULLNAME=$FILENAME
else
FULLNAME=alien:///alice/data/${YEAR}/${PERIOD}/000${RUN}/raw/${FILENAME}
fi
echo "FullName   ${FULLNAME}"
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
time aliroot -l <<EOI|tee execQA$RUN.log
EOF
.L $EXECFOLDER/ITSQArecoparam.C++
 ITSQArecoparam("${FULLNAME}",2,30); 
.q
EOI
time aliroot -l <<EOI|tee plot$RUN.log
.x $EXECFOLDER/PlotQASDD.C("ITS.QA.${RUN}.root");
.q
EOI
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
