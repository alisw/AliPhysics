#!/bin/bash

  [ -d qa ] || mkdir qa

  [ -z $1 ] && { echo "Usage: qa.sh <run_number>"; exit 1; }

  cd qa
  aliroot qa.C\($1\)
