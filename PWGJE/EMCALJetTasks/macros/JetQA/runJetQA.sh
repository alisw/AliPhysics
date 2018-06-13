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
#plotDir="$HOME/ALICE/JetQA/LHC15o/1872_highIR/TrainOutput"

#..........................................
#15o.......................................
#RUNLIST="246928 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246760 246759 246758 246757 246751 246750 246495 246493 246487 246434 246424 246271 246225 246222 246217 246115 246113 246089 246042 246037 246003 246001 245963 245954 245952 245949 245831 245829 245705 245702 245700 245683 246488 246087"
#plotDir="YourPlotDir"
#PREFIX="/alice/data/2015/LHC15o/000"
#SUFFIX="/pass1/PWGJE/Jets_EMC_PbPb/1872_20170310-0102/AnalysisResults.root"

#..........................................
#16j5.......................................
#RUNLIST="246991 246994"
#plotDir="YourPlotDir"
#PREFIX="/alice/"
#SUFFIX="YourTrainSuffix"

#..........................................
#15n.......................................
#RUNLIST="244340 244351 244355 244359 244364 244377 244411 244416 244418 244421 244453 244456 244480 244481 244482 244483 244484 244531 244540 244542 244617 244618 244619 244626 244627 244628"
#plotDir="YourPlotDir"
#PREFIX="/alice/data/2015/LHC15n/000"
#SUFFIX="YourTrainSuffix"

#..........................................
#17p.......................................
#RUNLIST="282343 282342 282341 282340 282314 282313 282312 282307 282306 282305 282304 282303 282302 282247 282230 282229 282227 282224 282206 282189 282147 282146 282126 282123 282122 282119 282118 282099 282098 282078 282051 282031 282030 282025"
#plotDir="YourPlotDir"
#PREFIX="/alice/data/2017/LHC17p/000"
#SUFFIX="YourTrainSuffix"

#..........................................
#17q.......................................
#RUNLIST="282367 282366 282365"
#plotDir="YourPlotDir"
#PREFIX="/alice/data/2017/LHC17q/000"
#SUFFIX="YourTrainSuffix"

downloadData=false
plotQA=false
generatePresentation=false

analysisFile="AnalysisResults.root"   # For pt-hard, set to "AnalysisResultsFinal.root"
referenceFile="AllRuns/AnalysisResults.root"

#-----------------------------------------
# Preparations
dir=$plotDir/"AllRuns"
if [ ! -d "$dir" ]; then
  mkdir -p $dir
  echo "created directory" $dir
fi
echo "output dir: $plotDir"

#-----------------------------------------
# Download Files Section
if [ "$downloadData" = true ]; then

  echo "-> Downloading data, per run..."

  for RUN in $RUNLIST
  do
    dir=$plotDir/$RUN
    if [ ! -d "$dir" ]; then
      mkdir $dir
      echo "created directory" $dir
    fi
    echo "- - - - - - - - - - - - - - - - -"
    echo "Downloading run " $RUN
    alien_cp alien://$PREFIX$RUN$SUFFIX $dir
  done
fi

#-----------------------------------------
# Plot QA of Files Section
if [ "$plotQA" = true ]; then

  echo "-> Executing plotting macro..."

  echo "- - - - - - - - - - - - - - - - -"
  echo "Plot merged reference run"
  python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/plotPWGJEQA.py -f $plotDir/$referenceFile -o $plotDir/"AllRuns"/QAoutput -i ".png"

  for RUN in $RUNLIST
  do
    echo "- - - - - - - - - - - - - - - - -"
    echo "Plotting run " $RUN
    python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/plotPWGJEQA.py -f $plotDir/$RUN/$analysisFile -o $plotDir/$RUN/QAoutput -r $plotDir/$referenceFile -i ".png"
  done
fi

#-----------------------------------------
# Generate powerpoint presentation Section
if [ "$generatePresentation" = true ]; then

  echo "-> Generating presentation..."
  python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/plotPowerpoint.py -r "$RUNLIST" -d "$plotDir"
fi

echo "Done!"
