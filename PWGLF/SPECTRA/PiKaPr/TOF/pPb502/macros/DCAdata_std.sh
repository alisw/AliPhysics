#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching DCAdata_std on $1 in background"
aliroot -b -q "DCAdata_std.C(\"$1\")" &> DCAdata.$1.log &
echo "DCAdata_std started"