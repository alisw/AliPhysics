#!/bin/bash

#useMB="$1"
#useExtended="$2"
#isMC="$3"

#if [ -z "$useExtended" ]; then
#  argument="$useMB"
#else
#  if [ -z "$isMC" ]; then
#    argument="$useMB,$useExtended"
#  else
#    argument="$useMB,$useExtended,$isMC"
#  fi
#fi
#
#. scripts/extract_signal.sh $argument
root -l -b -q Spectra.cc+
root -l -b -q Systematics.cc+
root -l -b -q Systematics.cc+'(true)'
root -l -b -q JoinSystematics.cc+
root -l -b -q Final.cc+
root -l -b -q BWFits.cc+
root -l -b -q BWFits.cc+'(true)'
root -l -b -q Final.cc+
root -l -b -q Ratio.cc+
