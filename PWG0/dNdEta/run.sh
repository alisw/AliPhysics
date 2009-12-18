#!/bin/bash

for I in CINT1A-ABCE-NOPF-ALL CINT1C-ABCE-NOPF-ALL CINT1-E-NOPF-ALL
#for I in CINT1B-ABCE-NOPF-ALL CINT1A-ABCE-NOPF-ALL CINT1C-ABCE-NOPF-ALL CINT1-E-NOPF-ALL
do
  root -b -q 'run.C(0, "/PWG0/jgrosseo/run", -1, 0, 0, 2, 2, "SAVE", "'$I'")'
done
  
