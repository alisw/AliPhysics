#!/bin/sh 
#############################################################################
# recpoints.sh. Front-end script to run reconstruction from the grid chunks
# Usage:
#    ./recpoints.sh <run_number> <pass_number>
#############################################################################
#
# modification history
# version 1.0  March 2010 
# author Yves Schutz CERN
#

# SET THE FOLLOWING PARAMETERS IF NEEDED: 
# ---------------------------------------
# ---------------------------------------
kill -9 `ps | grep aliroot | awk '{print $1}'`

export RUNNUM=$1
export PASS=$2
export DET=$3
PROGRAM=aliroot

[ -z $RUNNUM ] && { echo "Please provide a run number..."; exit 1; }
[ -z $PASS ] && { echo "Please provide a pass number..."; exit 1; }
[ -z $DET ] && { echo "Please provide a detector name..."; exit 1; }

[ ! -e "$HOME/.globus/usercert.pem" ] && { echo "FAILED: There is no certificate in $HOME/.globus"; exit 1; }

alien-token-init 
source /tmp/gclient_env_$UID;

[ ! "$?" -eq "0" ] && { echo "FAILED: Token creation failed"; exit 1; }

VERSION=1.0
TITLE="Standalone QA checking of Grid recpoints chunks. v$VERSION"

# Retrieve the list of chunks from AliEn.......
export BASEDIR="/alice/data/20*"
PATTERN="$RUNNUM/ESDs/pass${PASS}/*${RUNNUM}*/${DET}.RecPoints.root"
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

 
for filename in $CHUNKS; do
    filename="${filename//\"}"
    set -- "${filename%  /pass${PASS}*}" 
    set -- "${1##*/pass${PASS}/}"
    CHUNK="${1/\/${DET}*}"
    echo "Running QA for chunk $filename. Outputs will be stored in "$RUNNUM"/"$CHUNK
    [ -d $RUNNUM"/"$CHUNK ] && echo $RUNNUM"/"$CHUNK "exists" || mkdir -p $RUNNUM"/"$CHUNK
    cd $RUNNUM"/"$CHUNK
    [ -f ${DET}.QA.${RUNNUM}.root ] && rm ${DET}.QA.${RUNNUM}.root
    [ -f QAImageQA.${RUNNUM}.ps ] && rm QAImageQA.${RUNNUM}.ps
    [ -f QA.root ] && rm QA.root
    [ -f ${DET}.RecPoints.root ] && echo "file ${DET}.RecPoints.root already there" || alien_cp alien://$filename .
    filename="${filename/${DET}.RecPoints/galice}"
    [ -f galice.root ] && echo "file galice.root already there" || alien_cp alien://$filename .
    filename="${filename/galice/AliESDs}"
    [ -f AliESDs.root ] && echo "file AliESDs.root already there" || alien_cp alien://$filename .
$PROGRAM -b<<EOF
  .L $ALICE_ROOT/prod/cosmic/recpointqa.C++
  recpointqa()
  .q 
EOF
cd ..
done
#ls */Merged.QA.Data.root > merged.list
#outfile="Merged.QA.Data."$RUNNUM".root"
#$PROGRAM -b <<EOF
#.L $ALICE_ROOT/test/cosmic/MergeQAMerged.C
#MergeQAMerged("$outfile", "merged.list") ; 
#.q
#EOF
#rm -f merged.list
#$PROGRAM -b -q $ALICE_ROOT/test/cosmic/qasummary.C 
#$PROGRAM -b  $ALICE_ROOT/test/QA/menuQA.C
#cd ..
