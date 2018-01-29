#!/usr/bin/env bash

# create list of the figures in the local directory
# make set of the xxx.html pages
period="LHC16g1"

aliroot -b -q '$AliPhysics_SRC/PWGPP/TPC/macros/tpcMCValidation.C+("'$period'")'  | tee  tpcMCValidation.log

