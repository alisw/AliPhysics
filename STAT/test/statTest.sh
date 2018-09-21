#!/bin/bash

# statTest.sh -- test suite for STAT
#
# How to run all STAT tests:
#
#   alienv enter AliRoot/latest  # load AliRoot environment
#   cd <AliRoot_Build_Directory>
#   ctest --output-on-failure -R func_STAT_*
#
# Or extra verbose:
#
#   ctest --extra-verbose -R func_STAT_*
#
# Tests output will be printed only in case of failures.

source $ALICE_ROOT/libexec/alilog4bash.sh
if [[ ! $ALIROOT_SOURCE ]]; then
  echo "ALIROOT_SOURCE must be defined to the source directory of AliRoot."
  exit 1
fi

export ROOT_HIST=0

testAliTreePlayer() {
  cp $ALIROOT_SOURCE/STAT/test/AliTreePlayerTest.C .
  root -n -b -l <<\EOF 2>&1 | tee testAliTreePlayer.log
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
    alilog_error "statTest.testAliTreePlayer: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS
}

testAliTMinutiToolkitTestLinear() {
  cp $ALIROOT_SOURCE/STAT/test/AliTMinuitToolkitTestLinear.C .
  root -n -b -l <<\EOF 2>&1 | tee testAliTMinutiToolkitTestLinear.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./AliTMinuitToolkitTestLinear.C+(50,3)
EOF
  NPDF=$(ls -1 *.pdf | grep -c '\.pdf')
  if [[ $NPDF == 2 ]]; then
    alilog_success "statTest.testAliTMinutiToolkitTestLinear: All OK"
    exit 0
  else
    alilog_error "statTest.testAliTMinutiToolkitTestLinear: FAILED"
    exit 1
  fi
}


testAliDrawStyleTest() {
  cp $ALIROOT_SOURCE/STAT/test/AliDrawStyleTest.C .
  root -n -b -l <<\EOF 2>&1 | tee AliDrawStyleTest.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./AliDrawStyleTest.C+
EOF

  N_GOOD=$(grep -cE 'Ali.*OK' AliDrawStyleTest.log)
  N_BAD=$(grep -c "E-Ali" AliDrawStyleTest.log)

  TEST_STATUS=0
  if [[ $N_GOOD != 49 ]]; then
    alilog_error "statTest.AliDrawStyleTest: Test FAILED"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.AliDrawStyleTest: Test FAILED"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.AliDrawStyleTest: All OK"
  else
    alilog_error "statTest.: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS

}

testAliPainterTest() {
  cp $ALIROOT_SOURCE/STAT/test/AliPainterTest.C .
  root -n -b -l <<\EOF 2>&1 | tee AliPainterTest.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./AliPainterTest.C+
EOF

  N_GOOD=$(grep -cE 'Ali.*OK' AliPainterTest.log)
  N_BAD=$(grep -c "E-Ali" AliPainterTest.log)

  TEST_STATUS=0
  if [[ $N_GOOD != 5 ]]; then
    alilog_error "statTest.AliPainterTest: Test FAILED"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.AliPainterTest: Test FAILED"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.AliPainterTest: All OK"
  else
    alilog_error "statTest.: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS

}

testAliParserTest() {
  cp $ALIROOT_SOURCE/STAT/test/AliParserTest.C .
  root -n -b -l <<\EOF 2>&1 | tee AliParserTest.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./AliParserTest.C+
EOF

  N_GOOD=$(grep -cE 'Ali.*OK' AliParserTest.log)
  N_BAD=$(grep -c "E-Ali" AliParserTest.log)

  TEST_STATUS=0
  if [[ $N_GOOD != 30 ]]; then
    alilog_error "statTest.AliParserTest: Test FAILED"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.AliParserTest: Test FAILED"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.AliParserTest: All OK"
  else
    alilog_error "statTest.: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS

}

testTPCVolumeDemoTest() {
  cp $ALIROOT_SOURCE/STAT/Notebooks/TPCDataVolumeDemo.C .
  root -n -b -l <<\EOF 2>&1 | tee TPCVolumeDemo.C.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    .x ./TPCDataVolumeDemo.C
EOF

  N_GOOD=1
  N_BAD=0

  TEST_STATUS=0
  if [[ $N_GOOD != 1 ]]; then
    alilog_error "statTest.testTPCVolumeDemoTest: Test FAILED"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.testTPCVolumeDemoTest: Test FAILED"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.: All OK"
  else
    alilog_error "statTest.testTPCVolumeDemoTest: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS

}
testAliNDFunctionInterface() {
  cp $ALIROOT_SOURCE/STAT/Macros/AliNDFunctionInterface.cxx .
  cp $AliRoot_SRC/STAT/Macros/QAtrendingFitExample.C .
  root -n -b -l <<\EOF 2>&1 | tee AliNDFunctionInterface.cxx.log
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    AliDrawStyle::SetDefaults()
    AliDrawStyle::ApplyStyle("figTemplate");
    gStyle->SetOptTitle(1);
    .L AliNDFunctionInterface.cxx+
    .x QAtrendingFitExample.C
EOF
  N_GOOD=1
  N_BAD=0

  TEST_STATUS=0
  if [[ $N_GOOD != 1 ]]; then
    alilog_error "statTest.testAliNDFunctionInterface: Test FAILED"
    ((TEST_STATUS++))
  fi
  if [[ $N_BAD != 0 ]]; then
    alilog_error "statTest.testAliNDFunctionInterface: Test FAILED"
    ((TEST_STATUS+=2))
  fi
  if [[ $TEST_STATUS == 0 ]]; then
    alilog_success "statTest.: All OK"
  else
    alilog_error "statTest.testAliNDFunctionInterface: FAILED (code $TEST_STATUS)"
  fi
  exit $TEST_STATUS

}
[[ $1 ]] && $1
