#!/bin/bash

source /afs/cern.ch/alice/caf/caf-lxplus.sh -alien v4-18-12-AN

  [ -d qa_pp ] || mkdir qa_pp

  [ -z $1 ] && { echo "Usage: qa_pp.sh <run_number>"; exit 1; }

  cd qa_pp
 root.exe qa_pp.C\($1\)

  rm plot_macros/QAsym.proof.root
  ln -s /home/alishift/acrcaf/qa_pp/run$1.root plot_macros/QAsym.proof.root

