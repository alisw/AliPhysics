#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching TPCcalib on $1 in background"
aliroot -b -q "TPCcalib.C(\"$1\")" &> TPCcalib.$1.log &
echo "TPCcalib started"