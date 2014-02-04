#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching TOFmatchMC on $1 in background"
aliroot -b -q "TOFmatchMC.C(\"$1\")" &> TOFmatchMC.$1.log &
echo "TOFmatchMC started"