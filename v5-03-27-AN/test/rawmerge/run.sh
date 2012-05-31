#!/bin/bash

if ! test -e "sig/galice.root"; then
    mkdir -p sig
    cp makeSDigit.C Config.C sig
    cd sig
    aliroot -b -q makeSDigit.C'(117048,"/alice/data/2010/LHC10b/000117048/raw/10000117048041.30.root")'
    cd ..
fi
exit;
if ! test -e "bkg/galice.root"; then
    cp makeSDigit.C Config.C bkg
    mkdir -p bkg
    cd bkg
    aliroot -b -q makeSDigit.C'(117048,"/alice/data/2010/LHC10b/000117048/raw/10000117048041.40.root")'
    cd ..
fi

if ! test -e "merged/galice.root"; then
    mkdir merged
    cp -a sig/*.root merged
    cp merge.C rec.C merged
    cd merged
    aliroot -b -q merge.C
    aliroot -b -q rec.C
    cd ..
fi
