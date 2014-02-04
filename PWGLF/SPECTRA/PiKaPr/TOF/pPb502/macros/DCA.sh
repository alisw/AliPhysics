#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument needed"
    exit 1
fi
echo "launching DCAdata on $1 in background"
aliroot -b -q "DCA.C(\"$1\", 0)" &> DCA.$1.log &
echo "DCAdata started"