#! /bin/bash

if [ "$1" = "" ]; then
    echo "filename not specified"
    exit 1
fi

if [ "$2" = "" ]; then
    nevents=kMaxInt 
else
    nevents=$2
fi

echo "*********************************************"
echo " make libs"
echo "*********************************************"

aliroot -b -q MakeLibs.C &> makelibs.log

echo "*********************************************"
echo " starting TOFpid on $1"
echo " run over $nevents events"
echo "*********************************************"

for (( icharge = 0; icharge < 2; icharge++ )); do
    for (( ipart = 2; ipart < 5; ipart++ )); do
	echo aliroot -b -q "TOFpid.C(\"$1\", $ipart, $icharge, $ipart, kTRUE, kFALSE, kFALSE, -2., 2., $nevents)"
	aliroot -b -q "TOFpid.C(\"$1\", $ipart, $icharge, $ipart, kTRUE, kFALSE, kFALSE, -2., 2., $nevents)" &> "TOFpid.$ipart.$icharge.$ipart.$1.log" &

#	echo aliroot -b -q "TOFpid.C(\"$1\", $ipart, $icharge, $ipart, kFALSE, kFALSE, kFALSE, -2., 2., $nevents)"
#	aliroot -b -q "TOFpid.C(\"$1\", $ipart, $icharge, $ipart, kFALSE, kFALSE, kFALSE, -2., 2., $nevents)" &> "TOFpid.norapiditycut.$ipart.$icharge.$ipart.$1.log" &
	
    done
done

echo "*********************************************"
echo " TOFpid processes now running in background"
echo "*********************************************"

