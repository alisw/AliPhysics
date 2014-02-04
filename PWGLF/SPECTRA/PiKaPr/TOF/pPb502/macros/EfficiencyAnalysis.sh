#! /bin/bash

if [ "$2" == "" ]; then
    echo "two arguments needed"
    exit 1
fi

aliroot -b -q MakeLibs.C
./TOFmatchEff.sh $1
./TOFmatchEff.sh $2
./TOFmatchMC.sh $2
./TrackingEff.sh $2
