#!/bin/bash

# This scripts contains the functionality to download run-by-run AnalysisResults.root from AliEn, generate
# QA plots for each run, and plot the run-by-run QA into a powerpoint presentation.
#
# Set the below config variables as desired, then run this script from anywhere.
#
# Pre-reqs:
#    1) edit below config variables as necessary:
#       - RUNLIST: List of runs to iterate over.
#       - plotDir: Directory where data will be downloaded, and where QA will be output.
#       - downloadData: Flag specifying whether to download data.
#       - plotQA: Flag specifying whether to generate QA plots, to be located in $plotDir/RUNX/QAoutput
#       - generatePresentation: Flag specifying whether to generate powerpoint presentation, to be located in $plotDir
#       - analysisFile: Name of .root file to be analyzed for each run.
#       - referenceFile: Name of reference .root file (i.e. sum of all runs), for comparison to run-by-run plots (if plotQA enabled)
#       - PREFIX and SUFFIX: If downloading enabled, specify location on AliEn. The files should be located at PREFIX/RUN/SUFFIX.
#    2) load aliroot environment (and get alien token, if downloading data enabled)
#
# If option to download data is disabled (e.g. if pt-hard merging is done), the data files should be
# organized in the same file structure, i.e. a directory for each run, and $plotDir/TrainOutput should be the parent directory.
#
# Author: James Mulligan <james.mulligan@yale.edu>

RUNLIST="246928 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246760 246759 246758 246757 246751 246750 246495 246493 246487 246434 246424 246271 246225 246222 246217 246115 246113 246089 246042 246037 246003 246001 245963 245954 245952 245949 245831 245829 245705 245702 245700 245683 246488 246087"

plotDir="$HOME/ALICE/JetQA/LHC15o/1872_highIR/TrainOutput"

downloadData=false
plotQA=false
generatePresentation=false

analysisFile="AnalysisResults.root"   # For pt-hard, set to "AnalysisResultsFinal.root"
referenceFile="../AnalysisResultsReference.root"

PREFIX="/alice/data/2015/LHC15o/000"
SUFFIX="/pass1/PWGJE/Jets_EMC_PbPb/1872_20170310-0102/AnalysisResults.root"

if [ ! -d "$plotDir" ]; then
  mkdir "$plotDir"
fi
cd $plotDir
echo "output dir: $plotDir"

if [ "$downloadData" = true ]; then

  echo "Downloading data, per run..."

  for RUN in $RUNLIST
  do
    if [ ! -d "$RUN" ]; then
      mkdir "$RUN"
    fi
  done

  for RUN in $RUNLIST
  do
    alien_cp alien://$PREFIX$RUN$SUFFIX $RUN
  done

fi

if [ "$plotQA" = true ]; then

  echo "Executing plotting macro..."

  for RUN in $RUNLIST
  do
    python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/plotPWGJEQA.py -f $RUN/$analysisFile -o $RUN/QAoutput -r $referenceFile -i ".png"
  done

fi

if [ "$generatePresentation" = true ]; then

  echo "Generating presentation..."
  python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/plotPowerpoint.py -r "$RUNLIST" -d "$plotDir"

fi

echo "Done!"
