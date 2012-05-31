#!/bin/bash

RUN='123456'
PERIOD='LHC10b'
PASS='pass1'
YEAR='2010'
ALICE_ITS='$ALICE_ROOT/ITS'
TMPPLACE='/tmp'
TMPFOLDER='1'
EXECFOLDER='$HOME/macroQAshifter'
FILERUN=LHC10b.txt
FILEQATREND='FileQAtrend'
echo "file with the run numbers (Please insert the full path without environment variables)  :[${FILERUN}]"
read
if [ "$REPLY" != "" ]; then
export FILERUN="$REPLY"
echo "File    $FILERUN"
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
echo "local or lxplus (1=local 2=lxplus)   :[${TMPFOLDER}]"
read
if [ "$REPLY" != "" ]; then
export TMPFOLDER=$REPLY
fi
if [ "$TMPFOLDER" == "1" ]; then
export TMPPLACE='/tmp'
else
export TMPPLACE='/tmp/$USERNAME'
fi
echo "folder with macros     :[${EXECFOLDER}]"
read
if [ "$REPLY" != "" ]; then
export EXECFOLDER=$REPLY
echo "Folder:    $EXECFOLDER"
fi
time aliroot  -l <<EOI | tee trend.log 
.L $EXECFOLDER/TrendQASDD.C++
 TrendQASDD("${FILERUN}","FileQAtrend","${PASS}",$YEAR,"${PERIOD}") 
.q
EOI
gv SDDtrend$PERIOD$PASS.ps & 
echo "Done!!"
