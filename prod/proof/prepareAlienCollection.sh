#!/bin/sh
#############################################################################
# prepareCollection.sh. Script to prepare a raw-data chunks collection
# for a given run
# Usage:
#    ./prepareCollection.sh <run_number> <collection_file>
#############################################################################
#
# modification history
# version 1.0  2008/09/12 Cvetan Cheshkov

# SET THE FOLLOWING PARAMETERS IF NEEDED: 
# ---------------------------------------
YEAR=09
# ---------------------------------------

RUNNUM=$1

[ -z $RUNNUM ] && { echo "Please provide a run number..."; exit 1; }

OUTFILE=$2

[ -z $OUTFILE ] && { echo "Please provide an output filename..."; exit 1; }

[ ! -e "$HOME/.globus/usercert.pem" ] && { echo "FAILED: There is no certificate in $HOME/.globus"; exit 1; }

[ -e "/tmp/gclient_env_$UID" ] && { source /tmp/gclient_env_$UID; }
alien-token-init 

[ ! "$?" -eq "0" ] && { echo "FAILED: Token creation failed"; exit 1; }

# Retrieve the list of chunks from AliEn.......
BASEDIR="/alice/data/20"$YEAR
PATTERN="/raw/"$YEAR"0000"$RUNNUM"*0.root"
rm -f collection.tmp
gbbox find $BASEDIR $PATTERN | grep -v found | head --lines=-1 > collection.tmp

[ $(stat -c%s collection.tmp) -eq 0 ] && { echo "No chunks found for the given run"; exit 1; }

rm -f $OUTFILE
for ifile in `cat collection.tmp` ; do echo "alien://$ifile" >> $OUTFILE 2>&1; done
rm -f collection.tmp

echo `cat $OUTFILE | wc -l`" raw-data chunks are found and added to the collection file "$OUTFILE;

