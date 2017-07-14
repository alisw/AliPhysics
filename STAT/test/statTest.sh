#!/bin/bash

# statTest.sh -- test suite for STAT
#
# How to run all STAT tests:
#
#   alienv enter AliRoot/latest  # load AliRoot environment
#   cd <AliRoot_Build_Directory>
#   ctest --output-on-failure -R func_STAT_*
#   or extra verbose
#   ctest --extra-verbose -R func_STAT_*
# Tests output will be printed only in case of failures.

source $ALICE_ROOT/libexec/alilog4bash.sh
if [[ ! $ALIROOT_SOURCE ]]; then
  echo "ALIROOT_SOURCE must be defined to the source directory of AliRoot."
  exit 1
fi

testAliTreePlayer() {
  # MI comment: Keep log file for later log archiving and log parsing
  # TMPDIR=$(mktemp -d)
  # cd $TMPDIR
  cp $ALIROOT_SOURCE/STAT/test/AliTreePlayerTest.C .
  root -b -l <<\EOF &> testAliTreePlayer.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./AliTreePlayerTest.C+
EOF
  N_GOOD=$(grep -cE 'AliTreePlayerTest\..*Test.*OK' testAliTreePlayer.log)
  N_BAD=$(grep -c "E-" testAliTreePlayer.log)
  TEST_STATUS=0
  if [[ $N_GOOD != 4 ]]; then
    alilog_error "statTest.testAliTreePlayer: Invariant test failed"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.testAliTreePlayer: Invariant test failed"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.testAliTreePlayer: All OK"
  else
    alilog_success "statTest.testAliTreePlayer: FAILED (code $TEST_STATUS): full log follows"
    cat testAliTreePlayer.log
  fi
  # cd /
  # rm -rf $TMPDIR
  exit $TEST_STATUS
}

testAliTMinutiToolkitTestLinear() {
  #MI comment: Keep log file for later log archiving and log parsing
  #TMPDIR=$(mktemp -d)
  #cd $TMPDIR
  cp $ALIROOT_SOURCE/STAT/test/AliTMinuitToolkitTestLinear.C .
  root -b -l <<\EOF &> testAliTMinutiToolkitTestLinear
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    cout << gSystem->GetIncludePath() << "e che cz" << endl;
    .x ./AliTMinuitToolkitTestLinear.C+(50,3)
EOF
  NPDF=$(ls -1 *.pdf | grep -c '\.pdf')
  if [[ $NPDF == 2 ]]; then
    alilog_success "statTest.testAliTMinutiToolkitTestLinear: All OK"
    exit 0
  else
    alilog_error "statTest.testAliTMinutiToolkitTestLinear: FAILED: full log follows"
    cat testAliTMinutiToolkitTestLinear
    exit 1
  fi
  #cd /
  #rm -rf $TMPDIR
}

[[ $1 ]] && $1
