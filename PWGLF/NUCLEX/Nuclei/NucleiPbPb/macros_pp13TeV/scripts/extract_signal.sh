#!/bin/bash
useMB="$1"
useExtended="$2"
isMC="$3"

if [ -z "$useExtended" ]; then
  argument="$useMB"
else
  if [ -z "$isMC" ]; then
    argument="$useMB,$useExtended"
  else
    argument="$useMB,$useExtended,$isMC"
  fi
fi

rootver=$(root-config --version)
if (( ${rootver:0:1} < 6 )); then
  echo "ROOT version 6 is required for this analysis."
else
  root -b -l << EOF
.L src/RooGausExp.cxx+
.L src/RooGausDExp.cxx+
.L src/FitModules.cxx+
.L Signal.cc+g
Signal($argument)

EOF

fi
