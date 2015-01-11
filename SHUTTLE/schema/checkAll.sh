#!/bin/bash

export DETS="ACO AD0 EMC FMD HMP MCH MTR PHS CPV PMD SPD SDD SSD TOF TPC TRD T00 V00 ZDC GRP"

if [ -n "$1" ] 
then
	DETS=$*
fi

for I in $DETS
do
        echo $I
	aliroot -b -q TestDPs.C\(\"$I\",1418209800,1418212000\) > $I.out 2>&1
	echo -n "OK: "
	grep Retrieved $I.out | wc -l
	echo -n "Failed: "
	grep Cannot $I.out | wc -l
	echo
done

