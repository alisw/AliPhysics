#!/bin/bash
echo ===========================
echo PATH= $PATH 
echo ROOTSYS = $ROOTSYS
echo LD_LIBRARY_PATH= $LD_LIBRARY_PATH
echo ALICE_ROOT = $ALICE_ROOT
echo ALICE_PHYSICS = $ALICE_PHYSICS
echo ==========================
free
echo _____________________________________________
echo "HOME IS $HOME"
ls $HOME
length=`echo $HOME |wc -c`
if   (( $length >= 100 )) ;
then
     echo "WARNING: The home directory $HOME is longer than 100 char"
     OLDHOME=$HOME
     NEWHOME="/tmp/alien_home_dir.${ALIEN_PROC_ID}"
     echo "CHANGING HOME TO $NEWHOME"
     ln -s "$HOME" "$NEWHOME"
     export HOME=$NEWHOME
fi
echo _____________________________________________
echo ALICE_ROOT = $ALICE_ROOT

export PRODUCTION_METADATA="$ALIEN_JDL_LPMMETADATA"

###

#if [ "$1" = "OCDB" ]; then
#    echo "Setting env for generating OCDB.root"
#
#    export OCDB_SNAPSHOT_CREATE="kTRUE"
#    export OCDB_SNAPSHOT_FILENAME="OCDB.root"
#
#    touch OCDB.generating.job
#
#    shift
#fi

### check environment
echo ALICE_ROOT = $ALICE_ROOT

if [ "$ALIDPG_ROOT" = "" ]; then

    if [ -f "alidpg.tgz" ]; then

	echo "Using AliDPG from tarball"
	tar zxvf alidpg.tgz
	export ALIDPG_ROOT=`pwd`/AliDPG
	
    else
	
	echo "*! ERROR: ALIDPG_ROOT is not set!"
	echo "ERROR: ALIDPG_ROOT is not set!" > validation_error.message
	exit
	
    fi
fi

### dgpsim.sh
echo ALICE_ROOT = $ALICE_ROOT

DPGSIMSH=$ALIDPG_ROOT/MC/dpgsim.sh
if [ -f dpgsim.sh ]; then
    chmod +x dpgsim.sh
    DPGSIMSH=./dpgsim.sh
fi
echo ALICE_ROOT = $ALICE_ROOT
echo "Calling '$DPGSIMSH $*'"
$DPGSIMSH $*
error=$?
echo ALICE_ROOT = $ALICE_ROOT

if [ $error -ne 0 ]; then
    echo "*! Command '$DPGSIMSH $*' exited with error code $error"
#    echo "Command '$DPGSIMSH $*' exited with error code $error" > validation_error.message
    exit $error
fi

root -q -b Tagging.C
# One can run another scan with different PHOSTender settings
# Can not be done within one root session as PHOSTender modifies PHOS clusters in AOD
#root -q -b Tagging16s.C

### 

if [ ! -z "$NEWHOME" ]; then
    echo "DELETING $NEWHOME"
    export HOME=$OLDHOME
    rm -rf $NEWHOME
fi

exit 0
