#! /bin/bash

if [ "$2" == "" ]; then
    echo "two arguments needed"
    exit 1
fi

aliroot -b -q MakeLibs.C
