#!/bin/bash
#######################################################################
#
# This script runs the test for the AliAnalysisTaskFiltered class 
#
#   Macro to test functionality of the AnaliAnalysisTaskFiltered.
#   To be used within UnitTest suit
#   $ALICE_ROOT/../src/test/testAliAnalysisTaskFiltered/AliAnalysisTaskFilteredTest.sh
#   To test:
#   1.) CPU/Memory/Data volume
#   2.) Relative fracion of the information in exported trees
#   3.) Compression for points
#
#   Author of test:
#   marian.ivanov@cern.ch
###############################################################################
# Test environment has to be configured before
# To setup: input data path 
# Setup example config can be found in directory $ALICE_ROOT/../src/test/configuration
# E.g for GSI setup:
#
# export UnitTestConfig=$ALICE_ROOT/../src/test/configuration/configGSI.sh 
# $ALICE_ROOT/../src/test/testAliAnalysisTaskFiltered/AliAnalysisTaskFilteredTest.sh
#
###############################################################################
# Steps:
#  1. define vars
#  2. echo settings
#  3. runTask 
###############################################################################

source $UnitTestConfig

#get path to input list
    inputListfiles=$TestData_pPb
#get scale number for tracks
    filterT=${2-20}
#get scale number for V0s
    filterV=${3-2}
#get scale number of friends
    filterFriend=${4--1}
#get OCDB path (optional)
    OCDBpath=${5-"\"$OCDBPath_pPb\""}
#get max number of files 
    maxFiles=${6-"1000000"}
#get offset of first file
    offsetFile=${7-"0"}
#get max number of events
    maxEvents=${8-"30000000"}
#get offset of first event
    offsetEvent=${9-"0"}

# echo settings
    if [[ -f "$inputListfiles" ]] ; then
        inputListfiles="\"$inputListfiles\""
        echo "running with setup:"
        echo "ALICE_ROOT: $ALICE_ROOT"
        echo "inputListfiles: $inputListfiles"
        echo "scale tracks: $filterT"
        echo "scale V0s: $filterV"
        echo "scale Friends: $filterFriend"
        echo "OCDB path: $OCDBpath"
        echo "max Files: $maxFiles"
        echo "offset File: $offsetFile"
        echo "max Events: $maxEvents"
        echo "offset Event: $offsetEvent"
    else
        echo "inputListfiles not found: $inputListfiles"
        exit 1
    fi
#run FilterTask
    echo aliroot -l -b -q $ALICE_ROOT/../src/test/testAliAnalysisTaskFiltered/AliAnalysisTaskFilteredTest.C\($inputListfiles,$filterT,$filterV,$OCDBpath,$maxFiles,$offsetFile,$maxEvents,$offsetEvent\)
    aliroot -l -b -q $ALICE_ROOT/../src/test/testAliAnalysisTaskFiltered/AliAnalysisTaskFilteredTest.C\($inputListfiles,$filterT,$filterV,$filterFriend,$OCDBpath,$maxFiles,$offsetFile,$maxEvents,$offsetEvent\)
exit
