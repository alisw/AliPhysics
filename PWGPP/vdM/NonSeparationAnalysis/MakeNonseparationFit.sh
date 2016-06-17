#!/bin/bash

fill=$1;

if [ X$fill == X ]; then
    echo "USAGE: $0 fillNumber"
    exit 1
fi

#aliroot -b -q .x ExtractCTPScalers_${fill}.C
#aliroot -b -q .x MakeLumiRegion_${fill}.C

case $fill in
    4634)
	bcs="-1 344 464 827 1187 1558 1678 3177 3297"
	counter=0
	for bc in $bcs; do
##	    (aliroot -b -q .x MakeNonseparationFit_4634.C'(1,'$counter')' &&  aliroot -b -q .x MakePlots.C'(1,4634,1,'$bc')' ) &
	    (aliroot -b -q .x MakePlots.C'(1,4634,1,'$bc')' ) &
	    ((counter++))
	done
	wait
	;;
    *)
	echo fill $fill is not yet suppored
	;;
esac


#aliroot -b -q .x MakeNonseparationFit_${fill}.C'(kFALSE)'
#aliroot -b -q .x MakeNonseparationFit_${fill}.C'(kTRUE)'
