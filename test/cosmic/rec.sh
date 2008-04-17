#!/bin/sh
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
YEAR=08
DIALOG=${DIALOG=dialog}
GSHELL_ROOT=$ALIEN_ROOT/api
PATH=$PATH:$GSHELL_ROOT/bin
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
aliensh -c "gbbox find $BASEDIR $PATTERN" | head --lines=-1 > collection.tmp

[ $(stat -c%s collection.tmp) -eq 0 ] && { echo "No chunks found for the given run"; exit 1; }
list=`cat collection.tmp | awk '{printf("%s %s %s ",$1,"  .","0");}'`   # improve: put basename
totChunks=`cat collection.tmp | wc -l`; totChunks=`echo \($totChunks / 3\) | bc`
rm -f collection.tmp

tempfile=`tempfile 2>/dev/null` || tempfile=/tmp/test$$
trap "rm -f $tempfile" 0 1 2 5 15
$DIALOG --clear --no-cancel --title "$TITLE" \
        --ok-label OK --checklist "$totChunks chunks available. Select chunks for reconstruction" 18 74 10 \
        $list 2> $tempfile

CHUNKS=`cat $tempfile`" ."
echo "Selected chunks:"
echo $CHUNKS
echo

echo $CHUNKS | while read -d " " filename; do 
     filename=${filename//\"/} 
     CHUNK=`echo $filename | cut -d "/" -f 8 | cut -d "." -f 1,2 | cut -c 12-14,16-`

     echo "Running AliRoot reconstruction for chunk $filename. Outputs will be stored in "$RUNNUM"/"$CHUNK"."
     rm -rf   $RUNNUM"/"$CHUNK
     mkdir -p $RUNNUM"/"$CHUNK
     cd       $RUNNUM"/"$CHUNK
     aliroot -b -q ../../rec.C\(\"alien://$filename\"\) 2>&1 | tee rec.log
     cd ../..
done
