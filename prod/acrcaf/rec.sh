#!/bin/bash

[ -d reco/local ] || mkdir reco/local
[ -d reco/log   ] || mkdir reco/log

if [ "$1" == "-local" ] && [ ! -z $2 ]
then
    
  nev=${3:-10000}
  fev=${4:-0}

  cd reco/local
  mkdir -p run$2
  rm -f run$2/*
  cd run$2
  
  echo "Reconstructing into reco/local/run$2 and redirecting output to the file reco/local/run$2/stdout"

  unbuffer aliroot -q ../../rec.C\($2,$nev,$fev\) 2>&1 | tee stdout

else

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

  CAFROOT=v5-24-00b-caf
  BASE=
  
  if [ "$1" == "-ALICE_pro" ]
  then
    shift
 
    CAFROOT=ALICE_pro
    BASE=/afs/cern.ch/alice/library/afs_volumes/vol12
  fi

  if [ "$1" == "-ALICE_new" ]
  then
    shift
 
    CAFROOT=ALICE_new
    BASE=/afs/cern.ch/alice/library/afs_volumes/vol02
  fi

  if [ ! -z $BASE ]
  then
    echo "Setting ROOT and AliRoot to $CAFROOT"    

    export ROOTSYS=$BASE/root
    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

    export ALICE_ROOT=$BASE/AliRoot
    export ALICE_TARGET=`root-config --arch`
    export LD_LIBRARY_PATH=$ALICE_ROOT/lib/tgt_${ALICE_TARGET}:$LD_LIBRARY_PATH
    export PATH=$ALICE_ROOT/bin/tgt_${ALICE_TARGET}:$PATH
  fi 

  [ -z $1 ] && { echo "Usage: rec.sh [-local] <run_number>"; exit 1; }

  nev=${2:-10000}
  fev=${3:-0}

  cd reco
  aliroot -q run.C\($1,$nev,$fev,\"$CAFROOT\"\)

fi
