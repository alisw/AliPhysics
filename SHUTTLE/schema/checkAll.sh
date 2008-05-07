#!/bin/bash

export DETS="ACO EMC FMD HMP MCH MTR PHS CPV PMD SPD SDD SSD TOF TPC TRD T00 V00 ZDC GRP"

if [ -n "$1" ] 
then
	DETS=$*
fi

for I in $DETS
do
        echo $I
	aliroot -b -q TestDPs.C\(\"$I\",1210083077,1210083078\) > $I.out 2>&1
done

