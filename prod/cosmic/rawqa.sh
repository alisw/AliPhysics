#!/bin/sh 
#############################################################################
# rawqa.sh. Front-end script to run reconstruction from the grid chunks
# Usage:
#    ./rawqa.sh <run_number>
#############################################################################
#
# modification history
# version 1.0  July 2008 adapated from rec.C by M. Meoni
# author Yves Schutz CERN
#

# SET THE FOLLOWING PARAMETERS IF NEEDED: 
# ---------------------------------------
export YEAR=09
# ---------------------------------------

export RUNNUM=$1

[ -z $RUNNUM ] && { echo "Please provide a run number..."; exit 1; }

[ ! -e "$HOME/.globus/usercert.pem" ] && { echo "FAILED: There is no certificate in $HOME/.globus"; exit 1; }

#[ -e "/tmp/gclient_env_$UID" ] && { source /tmp/gclient_env_$UID; }
#echo 12==========================  $LD_LIBRARY_PATH
alien-token-init 
source /tmp/gclient_env_$UID;

[ ! "$?" -eq "0" ] && { echo "FAILED: Token creation failed"; exit 1; }

VERSION=1.0
TITLE="Standalone QA checking of Grid rawdata chunks. v$VERSION"

# Retrieve the list of chunks from AliEn.......
export BASEDIR="/alice/data/20"$YEAR
PATTERN="/raw/"$YEAR"0000"$RUNNUM"*0.root"
#aliensh -c "gbbox find $BASEDIR $PATTERN" | head --lines=-1 > collection.tmp
aliensh -c "gbbox find $BASEDIR $PATTERN" > collection.tmp

[ `ls -al collection.tmp | awk '{print $5}'` -eq 0 ] && { echo "No chunks found for the given run"; exit 1; }
rm -r collection.tmp2
#
for ifile in `cat collection.tmp | grep root` ; do printf $ifile" "\|" "0" " >> collection.tmp2 ; done
list=`cat collection.tmp2`
[ -e collection.tmp2 ] && { rm -f collection.tmp2 ; } 
totChunks=`cat collection.tmp | wc -l`
rm -f collection.tmp

tempfile=`tempfile 2>/dev/null` || tempfile=/tmp/test$$
trap "rm -f $tempfile" 0 1 2 5 15
dialog --clear --no-cancel --title "$TITLE" \
        --ok-label OK --checklist "$totChunks chunks available for run $RUNNUM (only the first 500 are shown). Select chunks for reconstruction" 18 80 10 \
        $list 2> $tempfile

CHUNKS=`cat $tempfile`
echo "Selected chunks:"
echo $CHUNKS
echo

#dialog --clear --no-cancel \
#        --ok-label OK --radiolist "Program to run:" 15 20 5 "aliroot -b" \| on alieve \| off 2> $tempfile
PROGRAM=aliroot #`cat $tempfile`

# 
for filename in $CHUNKS; do
     CHUNK=`basename $filename | cut -d "." -f 1,2`
     BEG=`expr index "$CHUNK" .`
     BEG=`expr $BEG - 4`
     SUBCHUNK=${CHUNK:$BEG}
     echo "Running QA for chunk $filename. Outputs will be stored in "$RUNNUM"/"$CHUNK".   $SUBCHUNK"
     [ -e $RUNNUM"/"$CHUNK ] && { rm -rf   $RUNNUM"/"$CHUNK ; }
     mkdir -p $RUNNUM"/"$CHUNK
     rm $RUNNUM"/"*.QA.$RUNNUM.$SUBCHUNK.root 
     rm $RUNNUM"/"QA.$SUBCHUNK.root 
     cd       $RUNNUM"/"$CHUNK
$PROGRAM -b <<EOF
.L $ALICE_ROOT/test/cosmic/rawqa.C+
rawqa($filename, $RUNNUM)
EOF

$PROGRAM -b <<EOF
AliQAManager * qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ; 
 qam.Merge(atoi(gSystem->Getenv("RUNNUM"))) ;
EOF
     rm *QA.$RUNNUM.root
     cd ..
done
ls */Merged.QA.Data.root > merged.list
outfile="Merged.QA.Data."$RUNNUM".root"
$PROGRAM -b <<EOF
.L $ALICE_ROOT/test/cosmic/MergeQAMerged.C
MergeQAMerged("$outfile", "merged.list") ; 
EOF
rm -f merged.list
#$PROGRAM -b -q $ALICE_ROOT/test/cosmic/qasummary.C 
#$PROGRAM -b  $ALICE_ROOT/test/QA/menuQA.C
cd ..
