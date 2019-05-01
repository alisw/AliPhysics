#!/bin/bash
#
# This script downloads Pt-hard bin data, and merges and scales appropriately using Pt-hard weights computed from histograms filled
# in AliAnalysisTaskPWGJEQA.
#
# The script has two modes, depending whether run-by-run QA is desired -- set flag doRunByRun.
#
# The script:
#   (1) Downloads output analysis files per Pt-hard bin and run (at the moment this assumes run-by-run output from train)
#   (2) Merges output files for each Pt-hard bin (per run, if run-by-run mode)
#   (3) Calls python script scalePtHardHistos.py to re-weight histograms in each Pt-hard bin file (per run, if run-by-run mode)
#   (4) Sums all weighted Pt-hard bins into a single final output file (per run, if run-by-run mode)
# To use:
#   Fill in appropriate fields in config. Then, execute script from anywhere,
#   and it will write out to folder TrainOutput in the specified outputDir.
#
# Author: James Mulligan <james.mulligan@yale.edu>

outputDir="$HOME/ALICE/JetQA/LHC16j5/1031_FullProduction"
doRunByRun=true

downloadFiles=true
mergeToPtHardBins=false
scaleHistograms=false
mergePtHardBins=false

useReferenceFile=true
referenceFile="$HOME/ALICE/JetQA/LHC16j5/1031_FullProduction/AnalysisResultsReference.root"

PtHardBins=20

# LHC15o
# RUNLIST="246810 246751 246994 245954 246991 246012 245952 245833 245949 246766 246765 246989 246003 246053 245705 246089 245831 246763 246276 245702 246185 246495 246052 246676 246225 246809 246153 246275 245700 246760 246493 246182 246049 246001 246675 245829 246759 246808 246222 245692 246152 246984 246758 246048 246807 246181 246087 246272 246757 246434 246151 246042 246488 246804 246948 246982 246846 246928 246115 246945 246180 246805 246217 246851 245683 246847 246750 246271 246487 246431 246844 246037 246178 246980 246845 246428 246424 245923 246113 246036 245064 246392 244982 244918 244975 244980 244983 245061 245066 245068 246390 246391 245152 245151 245146 245145 245232 245231 246648 246583 246575 246553 245793 245775 245738 245148 245411 245454 245259 245409 245507 245554 245452 245505 245407 245353 245450 245504 245545 245544 245446 245349 245501 245401 245543 245441 245347 245397 245542 245497 245346 245540 245439 245496 245396 245345 245535 245343 245785 245766 245759 245752 245731 245729 246871 246870 246867 246865 245410"

# LHC16j5 (EMCal good)
# RUNLIST="246945 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246760 246759 246758 246757 246751 246750 246495 246493 246488 246487 246434 246424 246272 246271 246225 246222 246217 246115 246113 246089 246087 246053 246052 246037 246003 246001 245954 245952 245949 245833 245831 245829 245705 245702 245700 245683"

data="sim"
year="2016"
period="LHC16j5"
trainName="Jets_EMC_pp_MC"
trainPWG="PWGJE"
trainNumber="1031_20170127-2234"
filename="AnalysisResults.root"

PREFIX="/alice/${data}/${year}/${period}/"
SUFFIX="/${trainPWG}/${trainName}/${trainNumber}/${filename}"

# Create outputDir and cd into it
if [ ! -d $outputDir/TrainOutput ]; then
  mkdir -p $outputDir/TrainOutput
fi
cd $outputDir/TrainOutput
echo "output dir: $outputDir/TrainOutput"

# (1) Download output analysis files per Pt-hard bin and run
if [ "$downloadFiles" = true ]; then
  # Create a folder for each Pt-hard bin and run number
  echo "Creating directory structure..."
  echo ""
  for bin in $(seq 1 $PtHardBins);
  do
    for RUN in $RUNLIST
    do
      if [ ! -d "$bin/$RUN" ]; then
        mkdir -p "$bin/$RUN"
      fi
    done
  done

  # Copy the files
  echo ""
  echo "Looking for file ${filename} in path ${PREFIX}PtHardBinXX/RunXXXXXX${SUFFIX}"
  echo "Downloading each AnalysisResults..."
  echo ""
  for bin in $(seq 1 $PtHardBins);
  do
    for RUN in $RUNLIST
    do
      echo "Copying $PREFIX$bin/$RUN$SUFFIX to TrainOutput/$bin/$RUN"
      echo "alien_cp alien://$PREFIX$bin/$RUN$SUFFIX $bin/$RUN"
      alien_cp alien://$PREFIX$bin/$RUN$SUFFIX $bin/$RUN
    done
  done
fi

# If run-by-run output not needed, proceed to merge runs, scale, and merge pt-hard bins
if [ "$doRunByRun" = false ]; then

  # (2) Merge output files for each Pt-hard bin
  if [ "$mergeToPtHardBins" = true ]; then
    for bin in $(seq 1 $PtHardBins);
    do
      cd $bin
      hadd AnalysisResultsPtHard$bin.root */AnalysisResults.root
      cd ..
    done
  fi

  # (3) Call python script to re-weight histograms in each Pt-hard bin file
  if [ "$scaleHistograms" = true ]; then
    if [ "$useReferenceFile" = true ]; then
      python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
    else
      python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py
    fi
  fi

  # (4) Sums all weighted Pt-hard bins into a single final output file.
  if [ "$mergePtHardBins" = true ]; then
    hadd AnalysisResultsFinal.root */AnalysisResultsPtHard*.root
  fi

else

  # Create dir structure for each run, and merge pt-hard bins for each run
  for bin in $(seq 1 $PtHardBins);
  do
    for RUN in $RUNLIST
    do
      if [ ! -d "$RUN/$bin" ]; then
        mkdir -p "$RUN/$bin"
      fi
      if [ -d "$bin" ]; then
        cp $bin/$RUN/AnalysisResults.root $RUN/$bin/AnalysisResultsPtHard$bin.root
      fi
    done
    if [ -d "$bin" ]; then
      rm -r $bin
    fi
  done

  # (3) Call python script to re-weight histograms in each Pt-hard bin file
  if [ "$scaleHistograms" = true ]; then
    for RUN in $RUNLIST
    do
      cd $RUN
      echo "Scaling Run $RUN"
      if [ "$useReferenceFile" = true ]; then
        python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
      else
        python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py
      fi
      cd ..
    done
  fi

  # (4) Sums all weighted Pt-hard bins into a single final output file, for each run.
  if [ "$mergePtHardBins" = true ]; then
    for RUN in $RUNLIST
    do
      cd $RUN
      echo "Merging Run $RUN"
      hadd AnalysisResultsFinal.root */AnalysisResultsPtHard*.root
      cd ..
    done
  fi

fi


