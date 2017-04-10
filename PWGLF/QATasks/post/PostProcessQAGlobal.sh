#!/bin/sh
####################################################
# Simple script to execute macros and produce output
# Author: Domenico Colella <domenico.colella@cern.ch>
# Date:   March  2017 
#
#  Requirement: the macro paths need to be defined
#  ex:  CODE=$ALICE_PHYSICS/PWGPP/analysisQA
#
#
#  sh process.sh
#
####################################################




CODEFOLDER="."                    # code folder path
outputfilename="LHC15o_IR"        # suffix to the pdf file name
presentationtitle="LHC15o{\_}IR"  # slide title [use latex text, if you want to use an underscore write {\_} (e.g. LHC15o{\_}IR)]
collsys="0"                       # 0 = PbPb, 1 = pp/pPb 
ismc="kFALSE"                     # kTRUE if MC production
log=AnalysisResults.log           # log file


#//////////////////////////
#  DEFINE USEFUL FUNCTIONS
#/////////////////////////
#__________________________________
# RUN POST-PROCESSING MACRO FOR PID
function pid(){
  echo " \n Processing PID outputs..." 2>&1 | tee -a $log
  echo "aliroot -l -q '$CODEFOLDER/MakeTrendingPIDQA.C\(\"AnalysisResults.root\",\"PIDqa\",\"all\",111111,\"PIDqaReport.pdf\"'" 2>&1 | tee -a $log
  aliroot -l -b -q $CODEFOLDER/MakeTrendingPIDQA.C\(\"AnalysisResults.root\", \"PIDqa\", \"all\", 111111, \"PIDqaReport.pdf\"\) 2>&1 | tee -a $log
  rm -rf ITSnSigmaPID.root TPC_TOFnSigmaPID.root TPCBasicnSigmaPID.root TPCV0nSigmaPID.root TOFnSigmaPID.root PIDqaReport.pdf trending.root
}
#_________________________________________
# RUN POST-PROCESSING MACRO FOR RESONANCES
function resonances(){
   echo " \n Processing resonance particles outputs..." 2>&1 | tee -a $log
}
#_______________________________________________________
# RUN POST-PROCESSING MACRO FOR SINGLE-STRANGE PARTICLES
function singlestrange(){
  echo " \n Processing single-strange particles outputs..." 2>&1 | tee -a $log
  echo " aliroot -l -q '$CODEFOLDER/PostProcessQAV0.C\(kTRUE,png)' " 2>&1 | tee -a $log
  aliroot -l -b -q $CODEFOLDER/PostProcessQAV0.C\(kTRUE,\"png\"\) 2>&1 | tee -a $log
}
#____________________________________________________
# RUN POST-PROCESSING MACRO FOR MULTI-STRANGE BARYONS
function multistrange(){
  echo " \n Processing multi-strange baryons outputs..." 2>&1 | tee -a $log
  echo " aliroot -l -q '$CODEFOLDER/PostProcessQAMultistrange.C\($collsys,$ismc,\"./\",png)' " 2>&1 | tee -a $log
  aliroot -l -b -q $CODEFOLDER/PostProcessQAMultistrange.C\($collsys,$ismc,\"./\",\"png\"\) 2>&1 | tee -a $log
}
#______________________________________
# MOVE ALL png FILES IN THE png/ FOLDER
function movepng(){
  if [ ! -d png ]; then mkdir png
  fi
  echo " \n Moving all the png files to png/ folder"
  count=`ls -1 *.png 2>/dev/null | wc -l`
  if [ $count != 0 ]
  then mv *".png" "./png/"
  else echo "Any .png file present!" 
  fi
}
#___________________________
# CREATE THE PDF USING LATEX
function presentation(){
  echo " \n STARTING TO PRODUCE THE PDF FILE"
  pdflatex "\providecommand{\ismc}{"$ismc"}\providecommand{\titlename}{"$presentationtitle"}\input{presentation.tex}"
  mv presentation.pdf "PWG-LF_QAanalysis_$outputfilename".pdf
}
#_________________
# CLEAN THE FOLDER
function cleanfolder(){
   rm presentation.aux presentation.log presentation.nav presentation.out presentation.snm presentation.toc
}

#//////////////////////////////////
#  CHOOSE FUNCTIONS TO BE EXECUTED
#//////////////////////////////////
#pid
#resonances
#singlestrange
multistrange
movepng
presentation
cleanfolder

















#_________________________
# DELETE SOURCE LATEX FILE
#rm presentation.tex


		 
