#! /bin/bash

aliroot -b - l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/FlatESDConverter.C++'("AliESDs.root",kTRUE,kTRUE)'
