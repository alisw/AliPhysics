#!/bin/bash
###########
# $Id$

export ALIHLT_USEPACKAGE=ROOT
export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/lib_$ALIHLT_USEPACKAGE

cd $ALIHLT_TOPDIR
make libs
