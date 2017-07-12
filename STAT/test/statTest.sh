#!/usr/bin/env bash
# Test tp be used in CTEST of make test of aliroot
#      see CMakeList.txt
# To execute test standalone:
#     ctest -R testAliTreePlayer --verbose
#     or
#     ( source $AliRoot_SRC/STEER/Utilities/alilog4bash.sh; source $AliRoot_SRC/STAT/test/statTest.sh; testAliTreePlayer )
#     ( source $AliRoot_SRC/STEER/Utilities/alilog4bash.sh; source $AliRoot_SRC/STAT/test/statTest.sh; testAliTMinutiToolkitTestLinear )
# assuming AliRoot_SRC variable is set

testAliTreePlayer(){
    #
    # see AliTreePlayerTest.C source code for the list of tests
    #     AliTreePlayerTest.log for details about test
    aliroot -b -q   $AliRoot_SRC/STAT/test/AliTreePlayerTest.C+ |tee AliTreePlayerTest.log
    cat AliTreePlayerTest.log | grep "AliTreePlayerTest\."
    nGOOD=`cat AliTreePlayerTest.log | grep "AliTreePlayerTest\." | grep  -c -i -e "TEST.*OK"`;
    nErr=`cat AliTreePlayerTest.log | grep -c "E-"`
    testStatus=0;
    if [ $nGOOD != 4 ] ; then
        alilog_error "statTest.testAliTreePlayer: Invariant test  failed. See log files AliTreePlayerTest.log"
        ((testStatus++))
    fi;
    if [ $nErr != 0 ] ; then
        alilog_error "statTest.testAliTreePlayer: Invariant test  failed. See log files AliTreePlayerTest.log"
        ((testStatus+=2))
    fi
    if [ $testStatus == 0 ] ; then
        alilog_success "statTest.testAliTreePlayer: All OK"
    else
        alilog_success "statTest.testAliTreePlayer: FAILED code %"
    fi;
    #exit 1; #
    exit $testStatus; # exist status tested in ctest
}

testAliTMinutiToolkitTestLinear(){
    aliroot -b -q  $AliRoot_SRC/STAT/test/AliTMinuitToolkitTestLinear.C+\(50,3\) |tee AliTMinuitToolkitTestLinear.log
    # for the moment test only that code did not fail
    # dead band on fit values to be checked in regression test (Elast search DB)?
    npdf=`ls *.pdf | grep -c pdf`
    if [ $npdf == 2 ] ;then
        alilog_success "statTest.testAliTMinutiToolkitTestLinear: All OK"
        exit 0;
    else
        alilog_error "statTest.testAliTMinutiToolkitTestLinear: FAILED"
        exit 1;
    fi
}