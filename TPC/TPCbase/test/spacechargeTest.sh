#!/bin/bash

# spaceChargeTest.sh -- test suite for AliTPCSpaceCharge3DDriftLine
#
# How to run all AliTPCSpaceCharge3DDriftLine tests:
#
#   alienv enter AliRoot/latest  # load AliRoot environment
#   cd <AliRoot_Build_Directory>

#
#   ctest --output-on-failure -R func_TPCbase_AliTPCSpaceCharge3DDriftLine
#  Or extra verbose:
#
#   ctest --extra-verbose -R func_TPCbase_AliTPCSpaceCharge3DDriftLine

#
# Tests output will be printed only in case of failures.

source $ALICE_ROOT/libexec/alilog4bash.sh
if [[ ! $ALIROOT_SOURCE ]]; then
  echo "ALIROOT_SOURCE must be defined to the source directory of AliRoot."
  exit 1
fi

export ROOT_HIST=0

testAliTPCSpaceCharge3DDriftLine() {
  cp $ALIROOT_SOURCE/TPC/TPCbase/test/AliTPCSpaceCharge3DDriftLineTest.C .
  root -n -b -l <<\EOF 2>&1 | tee testAliTPCSpaceCharge3DDriftLineTest.log
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  .x ./AliTPCSpaceCharge3DDriftLineTest.C+
EOF
  N_GOOD=$(grep -cE 'AliTPCSpaceCharge3DDriftLineTest.*Test.*OK.*' testAliTPCSpaceCharge3DDriftLineTest.log)
  N_BAD=$(grep -cE 'AliTPCSpaceCharge3DDriftLineTest.*Test.*FAILED.*' testAliTPCSpaceCharge3DDriftLineTest.log)
  TEST_STATUS=0
  if [[ $N_GOOD < 13 ]]; then
    # alilog_error "spacechargeTest.testAliTPCSpaceCharge3DDriftLine: Invariant test failed"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD > 0 ]]; then
    # alilog_error "spacechargeTest.testAliTPCSpaceCharge3DDriftLine: Invariant test failed"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "spacechargeTest.testAliTPCSpaceCharge3DDriftLine: All OK"
  else
    alilog_error "spacechargeTest.testAliTPCSpaceCharge3DDriftLine: FAILED (N_GOOD=$N_GOOD, N_BAD=$N_BAD)"
  fi

  exit $TEST_STATUS # test fix
 # exit 0  # This is hack - we need test to be running correctly or increase threshold  
}

[[ $1 ]] && $1
