#!/bin/bash

mkdir -p ppExport
coarse=1
if test $coarse -gt 0 ; then binning="coarse" ; else binning="full" ; fi

for sNN in 900 2760 7000 8000 ; do
    for trig in INEL V0AND ; do

	echo "======== sNN=${sNN} trig=${trig} ============="
	root -l -b -q Export.C\(\"pp\",${sNN},\"${trig}\",$coarse\)
	fname=`printf "pp_%04d_%s_%s_empirical.root" $sNN $trig $binning`
	if test -f $fname ; then
	    echo "Moving $fname ..."
	    mv plots/$fname ppExport/
	else
	    echo "$fname not found" > /dev/stderr
	    exit 1
	fi
    done
done 
    
