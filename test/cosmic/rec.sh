#!/bin/sh
#############################################################################
# rec.sh. Front-end script to run reconstruction from the grid chunks
# Usage:
#    ./rec.sh <run_number>
#############################################################################
#
# modification history
# version 1.0  2008/03/04 02:12  Marco Meoni

# SET THE FOLLOWING PARAMETERS IF NEEDED: 
# ---------------------------------------
YEAR=08
ALIENBIN=$ALIEN_ROOT/bin
DIALOG=${DIALOG=dialog}
PATH=$PATH:$ALIEN_ROOT/api/bin
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/alien/api/lib
# ---------------------------------------

RUNNUM=$1

if [ -z $RUNNUM ]; then
   echo "Please provide a run number..."
   exit
fi

[ -e $ALIENBIN/alien ] || { echo "Cannot find AliEn client: missing ALIEN_ROOT environment variable and/or AliEn installation..."; exit 1; }

VERSION=1.0
TITLE="Standalone reconstruction of Grid rawdata chunks. v$VERSION"
ALIEN="$ALIENBIN/alien login -exec"

# Retrieve the list of chunks from AliEn.......
BASEDIR="/alice/data/20"$YEAR
PATTERN=$YEAR"0000"$RUNNUM"*0.root"

$ALIEN " find $BASEDIR $PATTERN" > collection.tmp
[ $(stat -c%s collection.tmp) -eq 0 ] && { echo "No chunks found for the given run"; exit 1; }
list=`cat collection.tmp | awk '{printf("%s %s %s ",$1,"  .","0");}'`
rm -f collection.tmp

totChunks=`echo $list | wc -w`; totChunks=`echo \($totChunks / 3\) | bc`

tempfile=`tempfile 2>/dev/null` || tempfile=/tmp/test$$
trap "rm -f $tempfile" 0 1 2 5 15

$DIALOG --clear --no-cancel --title "$TITLE" \
        --ok-label OK --checklist "$totChunks chunks available. Select chunks for reconstruction" 18 74 10 \
        $list 2> $tempfile

CHUNKS=`cat $tempfile`" ."
echo "Selected chunks:"
echo $CHUNKS
echo

alien-token-init
. /tmp/gclient_env_$UID

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
