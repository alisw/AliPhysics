#!/bin/bash

cd "${ALICE_ROOT}"/../src/

# find TPC/ -maxdepth 1 -name '*.h' -or -name '*.cxx'

#Files=$( find TPC/ -name '*.C' -or -name '*.h' -or -name '*.cxx' )

# Files=(
#   'TPC/LandauTest.C'
#   #'TPC/Attic/AliTPCCalibTCF.h'
# )

# Files=(
#   'EVE/EveDet/AliEveEMCALSModuleData.h'
# )

Files=$( find TPC/ -maxdepth 1 -name '*.h' -or -name '*.cxx' )

while [[ $# -gt 0 ]] ; do
  case "$1" in
    -r) RestoreOnly=1 ;;
    -x) DoxygenOnly=1 ;;
    -d) Debug='--debug=debug' ;;
    -o) Stdout='-o' ;;
  esac
  shift 1
done

cd "${ALICE_ROOT}"/../src/doxygen

for F in ${Files[@]} ; do
  F="${ALICE_ROOT}/../src/${F}"
  [[ $DoxygenOnly != 1 ]] && git checkout "$F"
  [[ $RestoreOnly != 1 ]] && ./thtml2doxy.py $Stdout $Debug $( cat debug/includes.txt ) "$F"
done
