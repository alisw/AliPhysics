#!/bin/sh
####################################################
# Simple script to execute macros and produce output
# Author: Satyajit Jena <sjena@cern.ch>
# Date:   Mon Dec  9 13:28:09 CET 2013 
#
#  Requirement: the macro paths need to be defined
#  ex:  CODE=$ALICE_PHYSICS/PWGPP/analysisQA
#
#  Arguments for each macro should be change by
#  collecting proper information from respective
#  wagon owner.
#
#  sh process.sh <output.root> <aod-number>
#  ex: sh process.sh output.root 145 eps
#
####################################################

CODE="/Users/sjena/AnalysisQA/AnalysisQATrain/"
OUTPUT=$1       # output.root 
KEY=$2          # AOD number
log=running.log # log file
suffix=$3      # output type eps or png
INPUT=AnalysisResults.root


#_____________________________________________________________________
###running for CF Flow
echo " \n Processing Flow Outputs...." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processCFv2vsPt.C(\"$INPUT\",\"$suffix\",\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processCFv2vsPt.C\(\"$INPUT\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log


###running for Corr
# filename without .root at the end

echo " Processing Flow Outputs...."2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processMakeQA2pc.C(\"AnalysisResults\",\"$suffix\",\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processMakeQA2pc.C\(\"AnalysisResults\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log 

#_____________________________________________________________________
####running for DQ
echo " \n Processing J/Psi->ee Outputs...." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processJpsi2eeQAplots.C(\"jpsi_Default.root\",\"$suffix\",\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processJpsi2eeQAplots.C\(\"jpsi_Default.root\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log

#_____________________________________________________________________
###running for HFE Flow

echo " \n Processing Flow Outputs...." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processHFEQAtask.C(\"$INPUT\",\"$suffix\",\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processHFEQAtask.C\(\"$INPUT\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log

#_____________________________________________________________________
####running for JE
echo " \n Processing JE outputs... "2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processJETrackQA.C\(\"$INPUT\",\"$suffix\",10, 1, kFALSE, 0, \"$OUTPUT\")'" 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processJETrackQA.C\(\"$INPUT\",\"$suffix\",10,1,kFALSE,0,\"$OUTPUT\"\) 2>&1 | tee -a $log

####running for JE
echo " Processing JE outputs... "2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processJETriggerQA.C\(\"$INPUT\",\"$suffix\",0.2, 0.15, 0.3, 0,\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processJETriggerQA.C\(\"$INPUT\",\"$suffix\",0.2,0.15,0.3,0,\"$OUTPUT\"\) 2>&1 | tee -a $log


#_____________________________________________________________________
####running for LF
echo " \n Processing LF Multi Strange outputs..." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processMultistrangeQA.C\(0,0,\"./\",\"$INPUT\",\"$suffix\", \"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processMultistrangeQA.C\(0,0,\"./\",\"$INPUT\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log

#_____________________________________________________________________
#######running for GA
#
LM=Pi0IM_GammaTrackCorr_EMCAL_defaultCen0_100 #list name

echo " Processing GA  outputs..." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processDrawAnaCaloTrackQA.C\(\"$LM\",\"$INPUT\", \"$suffix\",0,\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processDrawAnaCaloTrackQA.C\(\"$LM\",\"$INPUT\",\"$suffix\",0,\"$OUTPUT\"\) 2>&1 | tee -a $log


####running for GA
#
LM=Pi0IM_GammaTrackCorr_EMCAL_EMCALCen0_100 #list name

echo " Processing GA  outputs..." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processDrawAnaCaloTrackQA.C\(\"$LM\",\"$INPUT\", \"$suffix\",0,\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processDrawAnaCaloTrackQA.C\(\"$LM\",\"$INPUT\",\"$suffix\",0,\"$OUTPUT\"\) 2>&1 | tee -a $log


####running for GA processProduceFastQA
CS=5080001022092970023220000000  # cut selection
EN=PbPb_2.76TeV
echo " Processing GA (processProduceFastQA) outputs..." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processProduceFastQA.C\(\"$INPUT\",\"$CS\",\"$suffix\",\"$EN\",\"\",\"AOD$KEY\",\"$OUTPUT\")'" 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processProduceFastQA.C\(\"$INPUT\",\"$CS\",\"$suffix\",\"$EN\",\"\",\"AOD$KEY\",\"$OUTPUT\"\) 2>&1 | tee -a $log 


#_____________________________________________________________________
###running for UD Flow

echo " \n Processing UD Outputs...." 2>&1 | tee -a $log
echo " aliroot -l -q '$CODE/processDrawUDQA.C(\"$INPUT\",\"$suffix\",\"$OUTPUT\")' " 2>&1 | tee -a $log
aliroot -l -b -q $CODE/processDrawUDQA.C\(\"$INPUT\",\"$suffix\",\"$OUTPUT\"\) 2>&1 | tee -a $log



		 
