#!/bin/bash

  if [ "$1" != "-force" ]
  then
    pgrep aliroot > /dev/null
    if [ "$?" -eq "0" ]
    then
      echo "Error: You can only start one rec.sh at a time. Please wait for the other session to terminate or close it."
      exit
    fi
  else
    shift
  fi

  [ -z $1 ] && { echo "Usage: rec.sh <run_number> <aaf_cluster> <num_events> <num_events_skip> <num_workers>"; exit 1; }

  PROOF="alice-caf.cern.ch"
  if [ ! -z $2 ]; then
    PROOF="$2"
  fi

  PROOFWK=""
  if [ ! -z $5 ]; then
    PROOFWK="$5"
  fi

  nev=${3:-1000}
  fev=${4:-0}

  cd reco
#   aliroot -q runProofLite.C\($1,$nev,$fev,\"VO_ALICE@AliRoot::$ALICE_LEVEL\"\)
  aliroot -q runProof.C\($1,$nev,$fev,\"aliprod@$PROOF\",\"$PROOFWK\"\)
