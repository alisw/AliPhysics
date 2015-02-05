#!/bin/bash

# Usage:
#   alidoxyimg.sh <macro>.C
#
# Generates <macro>.png image as the output of <macro>.
#
# Needs aliroot executable.

if [[ $1 == '' ]] ; then
  echo 'Usage:'
  echo '  alidoxyimg.sh <macro>.C'
  echo ''
  echo 'The output will be a file named <macro>.png.'
  exit 1
fi

if ! which aliroot > /dev/null 2>&1 ; then
  echo 'Error: aliroot executable not found in $PATH'
  exit 2
fi

base=${1%.*}
ext='png'

aliroot -b "$1" <<EOF
TVirtualPad *topPad = gPad;
while (topPad->GetMother() != topPad) {
  topPad = topPad->GetMother();
}
topPad->Print("$base.$ext");
EOF

if [[ ! -e "${base}.${ext}" ]] ; then
  echo "Error: output file ${base}.${ext} not generated!"
  exit 3
fi
