#!/bin/bash
###########
# $Id$

if test -z "$1"; then
 export ALIHLT_USEPACKAGE=ALIROOT
else
 export ALIHLT_USEPACKAGE=$1
fi
echo HLT ALIHLT_USEPACKAGE=$ALIHLT_USEPACKAGE

export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/lib_$ALIHLT_USEPACKAGE

#export ALIHLT_NOLOGGING=false
#export ALIHLT_DOMC=true
#export ALIHLT_ALIDETECT=true
