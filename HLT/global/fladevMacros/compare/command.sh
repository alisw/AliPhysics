#! /bin/bash

aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/CompareFlatESDs.C++'("flat.dat","normal.dat")' | tee compare.out
