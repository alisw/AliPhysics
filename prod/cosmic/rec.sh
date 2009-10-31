#!/bin/bash
#############################################################################
# rec.sh. Front-end script to run reconstruction from the grid chunks
# Usage:
#    ./rec.sh <run_number>
#############################################################################
#
# modification history
# version 1.0  2008/03/04 02:12  Marco Meoni
#
# version 1.1  2008/03/28 00:33  Marco Meoni. Used aliensh instead of alien client

# SET THE FOLLOWING PARAMETERS IF NEEDED: 
# ---------------------------------------
YEAR=09
DIALOG=`which dialog`
# ---------------------------------------

RUNNUM=$1

[ -z $RUNNUM ] && { echo "Please provide a run number..."; exit 1; }

[ ! -e "$HOME/.globus/usercert.pem" ] && { echo "FAILED: There is no certificate in $HOME/.globus"; exit 1; }

[ -e "/tmp/gclient_env_$UID" ] && { source /tmp/gclient_env_$UID; }
alien-token-init 

[ ! "$?" -eq "0" ] && { echo "FAILED: Token creation failed"; exit 1; }

VERSION=1.0
TITLE="Standalone reconstruction of Grid rawdata chunks. v$VERSION"

# Retrieve the list of chunks from AliEn.......
BASEDIR="/alice/data/20"$YEAR
PATTERN="/raw/"$YEAR"0000"$RUNNUM"*0.root"
gbbox find $BASEDIR $PATTERN | head -n 500 > collection.tmp

[ $(wc -l collection.tmp) -eq 0 ] && { echo "No chunks found for the given run"; exit 1; }
rm -f collection.tmp2
for ifile in `cat collection.tmp | head -n 500` ; do printf $ifile" "\|" "0" " >> collection.tmp2 ; done
list=`cat collection.tmp2`
rm -f collection.tmp2 
totChunks=`cat collection.tmp | wc -l`
rm -f collection.tmp

tempfile=`tempfile 2>/dev/null` || tempfile=/tmp/test$$
trap "rm -f $tempfile" 0 1 2 5 15
$DIALOG --clear --no-cancel --title "$TITLE" \
        --ok-label OK --checklist "$totChunks chunks available for run $RUNNUM (only the first 500 are shown). Select chunks for reconstruction" 18 80 10 \
        $list 2> $tempfile

CHUNKS=`cat $tempfile`
echo "Selected chunks:"
echo $CHUNKS
echo

$DIALOG --clear --no-cancel \
        --ok-label OK --radiolist "Program to run:" 15 20 5 "aliroot -b" \| on alieve \| off 2> $tempfile
PROGRAM=`cat $tempfile`
 
for filename in $CHUNKS; do
     filename=${filename//\"/}
     CHUNK=`basename $filename | cut -d "." -f 1,2`

     echo "Running AliRoot reconstruction for chunk $filename. Outputs will be stored in "$RUNNUM"/"$CHUNK"."
     rm -rf   $RUNNUM"/"$CHUNK
     mkdir -p $RUNNUM"/"$CHUNK
     cd       $RUNNUM"/"$CHUNK
     $PROGRAM -q $ALICE_ROOT/test/cosmic/rec.C\(\"alien://$filename\"\) 2>&1 | tee rec.log
     cd ../..
done
