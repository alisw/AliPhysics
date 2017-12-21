#!/bin/bash
rootver=$(root-config --version)
if (( ${rootver:0:1} < 6 )); then
  echo "ROOT version 6 is required for this analysis."
else
  root -b -l << EOF
.L src/RooGausExp.cxx+
.L src/RooGausDExp.cxx+
.L src/FitModules.cxx+
.x Signal.cc+g
EOF

fi

