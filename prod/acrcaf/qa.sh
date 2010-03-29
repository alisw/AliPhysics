#!/bin/bash

### ADDED THIS TEMPORARY TO FORCE A CERTAIN API SERVER ###
#@@@@@@@@     MAY NEED TO REMOVE IN FUTURE !!!!! @@@@@@@@
export alien_API_HOST=pcapiserv05.cern.ch
###########################################################

source /afs/cern.ch/alice/caf/caf-lxplus.sh -alien v4-18-Release.rec
#source /afs/cern.ch/alice/caf/caf-lxplus.sh -alien v4-19-04-AN

  [ -d qa ] || mkdir qa

  [ -z $1 ] && { echo "Usage: qa.sh <run_number>"; exit 1; }

  cd qa
  mkdir -p run$1
  rm -f run$1/*
  cp qa.C run$1/qa.C
  cd run$1
  echo "Running QA for run #$1 ..."
  echo "Output will be written to folder qa/run$1" 
  aliroot qa.C\($1\)
