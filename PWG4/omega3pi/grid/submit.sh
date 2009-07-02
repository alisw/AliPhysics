#!/bin/sh
#Usage: list.txt is the list of xml collections to process.

runs=`cat list.txt`

for run in $runs; do
    echo "submit Omega.jdl ${run}"
    submit Omega.jdl ${run}

done

