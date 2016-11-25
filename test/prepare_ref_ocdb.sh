#!/bin/bash -e
cd "$(dirname "$0")"
if [[ $OCDB_TEST_ROOT ]]; then
  find . -name "*.C" -exec perl -p -i -e 's|ALICE_ROOT/OCDB|OCDB_TEST_ROOT|' {} \;
else
  find . -name "*.C" -exec perl -p -i -e 's|OCDB_TEST_ROOT|ALICE_ROOT/OCDB|' {} \;
  OCDB_TEST_ROOT=$ALICE_ROOT/OCDB
fi
[[ -d $OCDB_TEST_ROOT/TPC ]] || { echo "Reference OCDB not found at $OCDB_TEST_ROOT! Export OCDB_TEST_ROOT to point to the correct OCDB path!"; false; }
