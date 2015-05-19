#!/bin/bash

Usage() {
    echo 'Usage: mprec2mille.sh [-s S=-200] <source1> <source2> ...'
    echo 'will convert source files or list of files <sourceI>'
    echo '<sourceI>.mille binaries'
    echo 'option -s S allows to control output chunk size in MB (S<0)'
    echo '            or number of tracks per chunk (S>0)            '
    echo 'S=0 : no splitting'
    exit 1
}

if [ $# -lt 1 ] ; then Usage ;fi

sz=-200

if [ "$1" == "-s" ] ; then
    shift 1 ;
    if [ $# -lt 2 ] ; then Usage ;fi
    sz="$1"
    shift 1 ;
fi

while [ $# -gt 0 ] ; do
inp=$1; 
shift 1 ;
#
 echo doing $inp
 aliroot -b -q MPRec2Mille.C+g\(\"${inp}\",0,${sz}\)
#
done

