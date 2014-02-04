#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching TOFmatchEff on $1 in background"
aliroot -b -q "TOFmatchEff.C(\"$1\")" &> TOFmatchEff.$1.log &
echo "TOFmatchEff started"