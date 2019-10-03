#! /bin/bash

if [ "$1" == "" ]; then
    echo "argument required"
    exit 1
fi

echo creating merged file: $1
ls TOFpid_*.root | tee list
aliroot -b -q "merger.C(\"list\", \"$1\", kFALSE)"



for ((icent = -1; icent < 7; icent++ )); do
    for (( icharge = 0; icharge < 2; icharge++ )); do
	for (( ipart = 2; ipart < 5; ipart++ )); do
	    echo "*** START FITTING ***"
	    echo "icent   = $icent"
	    echo "ipart   = $ipart"
	    echo "icharge = $icharge"
	    aliroot -b -q "defaultFit.C(\"$1\", $ipart, $icharge, $ipart, $icent)" &> defaultFit.$ipart.$icharge.$icent.$ipart.log &
	done
    done
    wait
done
	    
