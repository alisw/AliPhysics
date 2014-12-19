#!/usr/bin/env bash
#--------------------------------------------------------------------------------------------------
# Test suit for the utilities.sh 
#--------------------------------------------------------------------------------------------------
source ../alilog4bash.sh
source ../utilities.sh

test_guessRunData()
{
#
#  test functionalitu of guessRunNumber, guessYear, guessRunNumber from utilities.sh library
#  Few test use cases used. True values compared with "guessed" using utilities functionality
#  alilog_error generated in case of mismatch
#   guessRecoPass  to be added to the list ? 
    arrayPath=("/alice/data/2012/LHC12g/000188362/cpass1_tpc_validation/OCDB/spaceMap_4_4" \
               "alien:///alice/data/2012/LHC12a/000177011/cpass0_tpc_validation/OCDB/meanITSVertex.root" \
               "http://aliqatpc.web.cern.ch/aliqatpc/data/2010/LHC10e/pass4/000127729/TPC_dEdx_track_info.png"
        	"http://aliqatpc.web.cern.ch/aliqatpc/data/2010/LHC10e/pass4/meanTPCncl_vs_run.png" );
    yearTrue=("2012" "2012" "2010" "2010");
    periodTrue=("LHC12g" "LHC12a" "LHC10e" "LHC10e");
    runTrue=("188362" "177011" "127729" "");
    echo ===================================================================================================
    alilog "RUNNING TEST: test_guessRunData";
    echo ===================================================================================================
    #
    arrayLength=${#arrayPath[@]};	
    for (( i=0; i<${arrayLength}; i++ ));  do
        guessYear=`guessYear ${arrayPath[$i]}`
        guessPeriod=`guessPeriod ${arrayPath[$i]}`
        guessRun=`guessRunNumber ${arrayPath[$i]}`
        [[ "${guessYear}"   != "${yearTrue[$i]}"   ]] && alilog_error     "test_guessRunData: Year   ${guessYear}   ${yearTrue[$i]}"
        [[ "${guessYear}"   =  "${yearTrue[$i]}"   ]] && alilog_success   "test_guessRunData: Year   ${guessYear}   ${yearTrue[$i]}"
        [[ "${guessPeriod}" != "${periodTrue[$i]}" ]] && alilog_error     "test_guessRunData: Period   ${guessPeriod}   ${periodTrue[$i]}"
        [[ "${guessPeriod}" =  "${periodTrue[$i]}" ]] && alilog_success   "test_guessRunData: Period   ${guessPeriod}   ${periodTrue[$i]}"
        [[ "${guessRun}"    != "${runTrue[$i]}" ]] && alilog_error        "test_guessRunData: Run   ${guessRun}   ${runTrue[$i]}"
        [[ "${guessRun}"    =  "${runTrue[$i]}" ]] && alilog_success      "test_guessRunData: Run   ${guessRun}   ${runTrue[$i]}"
    done;
    echo ===================================================================================================
    alilog "END TEST: test_guessRunData";
    echo ===================================================================================================
}


