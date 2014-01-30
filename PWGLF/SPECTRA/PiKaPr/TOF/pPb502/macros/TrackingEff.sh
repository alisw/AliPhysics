#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching TrackingEff on $1 in background"
aliroot -b -q "TrackingEff.C(\"$1\")" &> TrackingEff.$1.log &
echo "TrackingEff started"