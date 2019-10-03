#!/bin/bash
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

if test -f EPOSLHC.tar.gz ; then
    tar -xzvf EPOSLHC.tar.gz
    (cd EPOSLHC && ./build.sh)
    if test $? -ne 0 ; then
	echo "Failed to build EPOS-LHC"
	exit 1
    fi
    cp -v EPOSLHC/libEPOSLHC.so .
    cp -v EPOSLHC/data/* .
fi 
export PRODUCTION_METADATA="$ALIEN_JDL_LPMMETADATA"

if test -f simrun.sh ; then 
    chmod 755 simrun.sh
    ./simrun.sh $@ 
    ret=$?
elif test -f Run.sh ; then 
    chmod 755 Run.sh
    ./Run.sh $@ 
    ret=$?
elif test -f run.sh ; then 
    chmod 755 Run.sh
    ./Run.sh $@ 
    ret=$?
else
    echo "Nothing to run in this directory!" 
    exit 1
fi
if test $ret -ne 0 ; then 
    echo "Job failed with exit status $ret"
    exit $ret
fi

exit 0
#
# EOF
#

    
