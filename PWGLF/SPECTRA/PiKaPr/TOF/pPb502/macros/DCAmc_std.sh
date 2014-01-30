#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching DCAmc_std on $1 in background"
aliroot -b -q "DCAmc_std.C(\"$1\")" &> DCAmc.$1.log &
echo "DCAmc_std started"