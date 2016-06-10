#!/bin/bash

fill=$1;

if [ X$fill == X ]; then
    echo "USAGE: $0 fillNumber"
    exit 1
fi

aliroot -b -q .x ExtractCTPScalers_${fill}.C
aliroot -b -q .x MakeLumiRegion_${fill}.C
aliroot -b -q .x MakeNonseparationFit_${fill}.C'(kFALSE)'
aliroot -b -q .x MakeNonseparationFit_${fill}.C'(kTRUE)'
