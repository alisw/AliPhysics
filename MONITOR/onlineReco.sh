#!/bin/bash

source ~/.bashrc

export BASE_DIR=/home/offline/onlineReco
export RECO_DIR=$BASE_DIR/reco

## enabling core dumps
ulimit -c 100000000000000

#----------------------------------------------------------------------
# check that no other instances of online reconstruction are running at the same time!

PIDFILE=$BASE_DIR/current.pid
if [ -f $PIDFILE ]
then
	OLDPID=`cat $PIDFILE`
	echo "Lock file $PIDFILE found. Checking process $OLDPID"
	kill -0 $OLDPID
	pidcheck=$?
	if [ "$pidcheck" -eq "1" ]
	then
		echo "Process $OLDPID done, removing file"
		rm -f $PIDFILE
	else
		echo "Process $OLDPID still running, exiting now!"
		exit 0;
	fi
fi

PROCNAME=`basename $0`
CURPID=`pgrep $PROCNAME`
echo "current process pid = " $CURPID 
echo $CURPID > $PIDFILE

cd $BASE_DIR

#------------------------------------------------------------------------
# init GRID environment

echo; echo 'Init GRID environment...'
#root_alien_setup

#------------------------------------------------------------------------
# Setting environment

echo;

export BUILD_DIR=$BASE_DIR/build
source $BUILD_DIR/SetAliRoot.sh

echo;

echo Root is: `root-config --version`
echo AliRoot is: `aliroot --version`
echo;

cd $BASE_DIR

if [ ! -e $RECO_DIR/log ]
then
	mkdir $RECO_DIR/log
fi

while [ 1 ] 
do
    cd $RECO_DIR/log
    aliroot -q $ALICE_ROOT/MONITOR/onlineReco.C\(\"listen\",\"$ALICE_ROOT/test/cosmic/rec.C\"\) | tee rec.log
#    mv rec.log log/run%%_rec.log
done

cd $BASE_DIR

# remove lock file
rm -f $PIDFILE

exit 0


