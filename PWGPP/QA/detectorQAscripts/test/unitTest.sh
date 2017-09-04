#!/usr/bin/env bash
#--------------------------------------------------------------------------------------------------
# Unit test for the detector QA script
# Currently implemented for the T0

#--------------------------------------------------------------------------------------------------
# Example to  run the test: 
#
# ( source $AliPhysics_SRC/PWGPP/QA/detectorQAscripts/test/unitTest.sh;  setupUnitTestMI; testQAT0; )
#
source $ALICE_ROOT/libexec/alilog4bash.sh
source $ALICE_ROOT/libexec/utilities.sh
#
# Configuration part:
#
setupUnitTestMI(){
    #
    #
    #
    export unitTestDir="${NOTES}/QA/ATO-102/data/UnitTest"
    export configFile="${NOTES}/QA/ATO-102/code/runQA-submitting.config"
    echo /lustre/nyx/alice/alienQA/alice/data/2013/LHC13e/000196309/pass2/QA_merge_archive.zip  > input.list
}

testQAT0(){
#
# UnitTest T0 
# To check success number of the output gif files used

    $ALICE_PHYSICS/PWGPP/QA/scripts/runQA.sh inputList="input.list" includeDetectors="T0" configFile="${configFile}"  workingDirectory="${unitTestDir}/T0/" outputDir="${unitTestDir}/T0/" logDirectory="${unitTestDir}/T0/"
    ngifT0=`find $wdir -iname "*gif" | grep -c gif`
    [ $ngifT0 = 50 ] && export unitTestT0=OK;
    [ $ngifT0 != 50 ] && export unitTestT0=FAILED;
    [ $ngifT0 = 50 ] && alilog_success      "unitTest.sh testQAT0 OK"
    [ $ngifT0 != 50 ] && alilog_error      "unitTest.sh testQAT0 FAILED"
}





