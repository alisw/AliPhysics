#!/bin/bash
###########
# $Id$

for i in ROOT ALIROOT STANDALONE; do
 export ALIHLT_USEPACKAGE=$i
 export ALIHLT_LIBDIR=$ALIHLT_TOPDIR/lib_$ALIHLT_USEPACKAGE
 echo  $ALIHLT_LIBDIR
 cd $ALIHLT_TOPDIR
 #make libs
 cd ..
done

