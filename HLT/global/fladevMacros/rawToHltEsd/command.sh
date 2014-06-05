#! /bin/bash

aliroot -b -q -l $ALICE_ROOT/HLT/exa/recraw-local.C'("raw.root","local://$OCDB10",-1,-1,"HLT","chains=GLOBAL-esd-converter ignore-hltout")' 2>&1|tee recraw-local.log

