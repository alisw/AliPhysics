#!/bin/bash
###########
# $Id$

for i in ROOT ALIROOT STANDALONE; do
 export ALIHLT_USEPACKAGE=$i
 export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/lib_$USEPACKAGE

 cd $ALIHLT_TOPDIR
 make clean
 cd ..
done
