#!/bin/sh
# $Id$
#
# by I. Hrivnacova, IPN Orsay

# Script to replace transformation/svmap data files 
# with generated ones.
# The old files are kept with *.old extension
#
# Usage: reset_data which_data
#                   which_data = transform, svmap 

if [ $# != 1 ] ; then
    echo " Usage: reset_data which_data"
    echo "                   which_data = transform, svmap"
    exit 1
fi

if [ "$1" != "transform" -a "$1" != "svmap" ] ; then
    echo " Usage: reset_data which_data"
    echo "                   which_data = transform, svmap"
    exit 1
fi

DATATYPE=$1

for builder in st1 st1V2 st2 slat trigger
do
  if [ -f $DATATYPE"_"$builder".dat.out" ] ; then
    # backup old file
    if [ -f $DATATYPE"_"$builder".dat" ] ; then
      mv $DATATYPE"_"$builder".dat" $DATATYPE"_"$builder".dat.old"
    fi  
    # activate new file
    mv $DATATYPE"_"$builder".dat.out" $DATATYPE"_"$builder".dat"
  echo ".. reset "$DATATYPE"_"$builder".dat" 
  fi
done  
