#! /bin/bash

if [ "$2" == "" ]; then
    echo "two arguments needed"
    exit 1
fi

aliroot -b -q MakeLibs.C
./DCAdata_std.sh $1
./DCAmc_std.sh $2
