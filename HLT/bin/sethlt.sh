#!/bin/bash
###########
# $Id$

export CVS_RSH=ssh
export CVSROOT=$USER@kjekspc1.fi.uib.no:/cvs/hltcvs
#export CVSIGNORE="lib src hough comp exa programs misc trigger sim 
#AliHLTCompCint.h AliHLTCompCint.cxx AliHLTHoughCint.h 
#AliHLTHoughCint.cxx AliHLTCint.cxx AliHLTCint.h lib_ROOT lib_ALIROOT
#AliHLTMiscCint.cxx AliHLTMiscCint.h"

export ALIHLT_USEPACKAGE=ALIROOT
#export ALIHLT_USEPACKAGE=ROOT
#export ALIHLT_USEPACKAGE=STANDALONE

export ALIHLT_BASEDIR=$HOME/work/hlt
export ALIHLT_TOPDIR=$ALIHLT_BASEDIR/level3code
export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/lib_$ALIHLT_USEPACKAGE

export ALIHLT_NOLOGGING=false
export ALIHLT_DOMC=true
export ALIHLT_ALIDETECT=true
export ALIHLT_ROWHOUGH=false
export ALIHLT_MLUCDIR=/usr/local/kip/MLUC

#export ALIHLT_DATADIR=/mnt/local/alidata/head
#export ALIHLT_TRANSFORMFILE=$ALIHLT_DATADIR/l3transform.config
#export ALIHLT_GEOPATH=$ALIDATADIR

if test -z "$LD_LIBRARY_PATH"; then
  export LD_LIBRARY_PATH=$ALIHLT_MLUCDIR/lib:$ALIHLT_LIBDIR
elif test -z "`echo $LD_LIBRARY_PATH | grep $ALIHLT_MLUCDIR/lib`"; 
then
  export LD_LIBRARY_PATH=$ALIHLT_MLUCDIR/lib:$ALIHLT_LIBDIR:$LD_LIBRARY_PATH 
fi
