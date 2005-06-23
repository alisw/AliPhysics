#!/bin/sh
# $Id$
#
# by I. Hrivnacova, IPN Orsay

# Script to compare existing transformation/svmap data files 
# with generated ones.
#
# Usage: compare_data which_data
#                     which_data = transform, svmap 

if [ $# != 1 ] ; then
    echo " Usage: compare_data which_data"
    echo "                     which_data = transform, svmap"
    exit 1
fi

if [ "$1" != "transform" -a "$1" != "svmap" ] ; then
    echo " Usage: reset_data which_data"
    echo "                   which_data = transform, svmap"
    exit 1
fi

DATATYPE=$1

for builder in st1 st1V2 st2 st2V2 slat trigger
do
  if [ -f $DATATYPE"_"$builder".dat.out" ] ; then
    echo "Comparing $DATATYPE"_"$builder".dat $DATATYPE"_"$builder.dat.out ...
    diff $DATATYPE"_"$builder.dat $DATATYPE"_"$builder.dat.out
  fi
done  
