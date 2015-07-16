#!/usr/bin/env bash
#--------------------------------------------------------------------------------------------------
# Test suit for the utilities.sh 
#--------------------------------------------------------------------------------------------------
# To run the test: 

# ( source $ALICE_PHYSICS/../src/PWGPP/scripts/test/utilities.test.sh; test_guessRunData)
# ( source $ALICE_PHYSICS/../src/PWGPP/scripts/test/utilities.test.sh; test_OCDBYear)
# ( source $ALICE_PHYSICS/../src/PWGPP/scripts/test/utilities.test.sh; test_paranoidCopy )

source $ALICE_PHYSICS/../src/PWGPP/scripts/alilog4bash.sh
source $ALICE_PHYSICS/../src/PWGPP/scripts/utilities.sh

test_guessRunData()
{
#
#  test functionalitu of guessRunNumber, guessYear, guessRunNumber from utilities.sh library
#  Few test use cases used. True values compared with "guessed" using utilities functionality
#  alilog_error generated in case of mismatch
#   guessRecoPass  to be added to the list ? 
    arrayPath=("/alice/data/2012/LHC12g/000188362/cpass1_tpc_validation/OCDB/spaceMap_4_4" \
               "alien:///alice/data/2012/LHC12a/000177011/cpass0_tpc_validation/OCDB/meanITSVertex.root" \
               "http://aliqatpc.web.cern.ch/aliqatpc/data/2010/LHC10e/pass4/000127729/TPC_dEdx_track_info.png" \
                "/hera/alice/local/filtered/toMerge/alice/data/2010/LHC10b/000117086/pass4/AOD/017/FilterEvents_Trees.root" \
        	"http://aliqatpc.web.cern.ch/aliqatpc/data/2010/LHC10e/pass4/meanTPCncl_vs_run.png" \
          "/hera/alice/local/filtered/alice/data/2010/LHC10e/000127712/ESDs/pass2/140/root_archive.zip"
               );
    yearTrue=("2012" "2012" "2010" "2010" "2010" "2010");
    periodTrue=("LHC12g" "LHC12a" "LHC10e"  "LHC10b" "LHC10e" "LHC10e" );
    runTrue=("188362" "177011" "127729"  "117086" "" "127712");
    passTrue=("cpass1_tpc_validation" "cpass0_tpc_validation" "pass4" "pass4_AOD" "pass4" "pass2_lego140");
    echo ===================================================================================================
    alilog "RUNNING TEST: test_guessRunData";
    echo ===================================================================================================
    #
    arrayLength=${#arrayPath[@]};	
    for (( i=0; i<${arrayLength}; i++ ));  do
	guessRunData ${arrayPath[$i]}
        guessYear=$year
        guessPeriod=$period
        guessRun=$runNumber
        guessPass=$pass
        [[ "${guessYear}"   != "${yearTrue[$i]}"   ]] && alilog_error     "test_guessRunData: Year   ${guessYear}   ${yearTrue[$i]} ${arrayPath[$i]}"
        [[ "${guessYear}"   =  "${yearTrue[$i]}"   ]] && alilog_success   "test_guessRunData: Year   ${guessYear}   ${yearTrue[$i]}"
        [[ "${guessPeriod}" != "${periodTrue[$i]}" ]] && alilog_error     "test_guessRunData: Period   ${guessPeriod}   ${periodTrue[$i]}"
        [[ "${guessPeriod}" =  "${periodTrue[$i]}" ]] && alilog_success   "test_guessRunData: Period   ${guessPeriod}   ${periodTrue[$i]} ${arrayPath[$i]}"
        [[ "${guessRun}"    != "${runTrue[$i]}" ]] && alilog_error        "test_guessRunData: Run   ${guessRun}   ${runTrue[$i]}   ${periodTrue[$i]} ${arrayPath[$i]}"
        [[ "${guessRun}"    =  "${runTrue[$i]}" ]] && alilog_success      "test_guessRunData: Run   ${guessRun}   ${runTrue[$i]}"
        [[ "${guessPass}"    != "${passTrue[$i]}" ]] && alilog_error        "test_guessPassData: Pass   ${guessPass}   ${passTrue[$i]}    ${periodTrue[$i]} ${arrayPath[$i]}"
        [[ "${guessPass}"    =  "${passTrue[$i]}" ]] && alilog_success      "test_guessPassData: Pass   ${guessPass}   ${passTrue[$i]}"
    done;
    echo ===================================================================================================
    alilog "END TEST: test_guessRunData";
    echo ===================================================================================================
}

test_OCDBYear(){
    #
    #  we had recently problem with this function in the StandardQA
    #  some slashes disapeared
    #  !!!!  In the guessRunData ocdbStorage is overwritten  !!!  
    #
    ocdbStorage="local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB/"
    newYear=2012
    newocdbStorage=$(setYear ${newYear} ${ocdbStorage})    
    newYear=2010
    newocdbStorage=$(setYear ${newYear} ${newocdbStorage})    
    [[  "${ocdbStorage}"  =   "${newocdbStorage}" ]] && alilog_success      "testOCDBYear: Pass"
    [[  "${ocdbStorage}"  !=   "${newocdbStorage}" ]] && alilog_error      "testOCDBYear: FAILED $ocdbStorage => $newocdbStorage"
}


test_paranoidCopy(){
    #
    # 1.) test the bug deleting of directoris
    #
    wdir=`pwd`;
    mkdir -p $wdir/testdirParanoidCopySource
    mkdir -p $wdir/testdirParanoidCopyDest
    echo test1 > $wdir/testdirParanoidCopySource/test1.txt
    echo test2 > $wdir/testdirParanoidCopySource/test2.txt
    #
    paranoidCp $wdir/testdirParanoidCopySource $wdir/testdirParanoidCopyDest > /dev/null
    paranoidCp $wdir/testdirParanoidCopySource $wdir/testdirParanoidCopySource > /dev/null
    #
    resultDiff=`diff -r  testdirParanoidCopySource testdirParanoidCopyDest/testdirParanoidCopySource`
    hrlen=${#resultDiff}
    [ $hrlen = 0 ]   && alilog_success     "test_paranoidCopy"	
    [ $hrlen != 0 ] && alilog_error      "test_paranoidCopy"	
}
