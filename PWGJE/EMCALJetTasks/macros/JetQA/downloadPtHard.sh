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
#         Eliane Epple   <eliane.epple@yale.edu>

outputDir="/YourOutPutDir"
###########  DEFAULT SETTINGS - Never Comment out ###########
doRunByRun=false
downloadFiles=false
mergeToPtHardBins=false
scaleHistograms=false
mergePtHardBins=false
useReferenceFile=false

#. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#In case you want to split your statistic into 2 sub-samples
mergeToPtHardBinsSubSample=false
scaleHistogramsSubSample=false
mergePtHardBinsSubSample=false

###########  DOWNLOAD THE FILES  ###########
#downloadFiles=true

###########  FOR REFERENCE FILE  ###########
##merge all runs into the 20 pT hard bins
#mergeToPtHardBins=true
#. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#mergeToPtHardBinsSubSample=true   #For half of the statistic

#--> AT THIS POINT WOULD BE GOOD TO BACKUP YOUR DATA <--
##scale the merged pT hard bins and then sum them
#useReferenceFile=true
#scaleHistograms=true
#mergePtHardBins=true
#. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#In case you want to split your statistic into 2 sub-samples
#scaleHistogramsSubSample=true
#mergePtHardBinsSubSample=true

###########  FOR SINGLE RUN FILES  ###########
##scale run-by-run the pT hard bins and then sum them - within one run
#doRunByRun=true
#scaleHistograms=true
#mergePtHardBins=true
#useReferenceFile=true

#----------------------------------------------------------------------------
#Period info
year="2018"
period="LHC18b8_fast"
declare -a RUNLIST=(282025 282031 282051 282078 282098 282099 282118 282119 282122 282123 282126 282146 282147 282189 282206 282224 282227 282229 282230 282247 282302 282303 282304 282305 282306 282307 282312 282313 282314 282340 282341 282342 282343 282365 282366 282367)

#year="2016"
#period="FILTER_PbPb_200_LHC16j5a"
#declare -a RUNLIST=(246945 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246760 246759 246758 246757 246751 246750 246495 246493 246488 246487 246434 246424 246272 246271 246225 246222 246217 246115 246113 246089 246087 246053 246052 246037 246003 246001 245954 245952 245949 245833 245831 245829 245705 245702 245700 245683)

NoOfRuns=${#RUNLIST[@]}
NoOfRunsHalf=$(expr $NoOfRuns / 2)
echo "There are ${NoOfRuns} number of runs"
echo "There are ${NoOfRunsHalf} half number of runs"
#----------------------------------------------------------------------------
#Train info
PtHardBins=20
data="sim"
trainName="Jets_EMC_pp_MC"
trainPWG="PWGJE"
trainNo="1581"
trainNoTime="_20180815-2319"   #No 1581
filename="AnalysisResults.root"

PREFIX="/alice/${data}/${year}/${period}/"
SUFFIX="/${trainPWG}/${trainName}/${trainNo}${trainNoTime}/${filename}"
referenceFile="$outputDir/ScaleFactors_Train1574.root"

#----------------------------------------------------------------------------
# Create outputDir
folderName="Train_"$trainNo

echo "output dir: $outputDir/$folderName"
if [ "$useReferenceFile" = true ]; then
  echo "Reference file: $referenceFile"
fi

if [ ! -d "${outputDir}/${folderName}" ]; then
  echo "Creating main directory..."
  mkdir "${outputDir}/${folderName}"
fi

# (1) Download output analysis files per Pt-hard bin and run
if [ "$downloadFiles" = true ]; then
  # Create a folder for each Pt-hard bin and run number

  echo "Creating directory structure..."
  echo ""
  for bin in $(seq 1 $PtHardBins);
  do
    for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
    do
      if [ ! -d "${outputDir}/${folderName}/$bin/${RUNLIST[$RUN]}" ]; then
        echo "Creating Folder ${outputDir}/${folderName}/$bin/${RUNLIST[$RUN]}"
        mkdir -p "${outputDir}/${folderName}/$bin/${RUNLIST[$RUN]}"
      fi
    done
  done

  # Copy the files
  echo ""
  echo "Looking for file ${filename} in path ${PREFIX}PtHardBinXX/RunXXXXXX${SUFFIX}"
  echo "Downloading each AnalysisResults..."
  echo ""
  for bin in $(seq 1 $PtHardBins);
  #for bin in $(seq 10 20);
  do
    echo "-------------------------"
    echo "pT hard bin No "$bin
    for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
    do
      echo "......"
      echo "Copying $PREFIX$bin/$RUN$SUFFIX to ...../$folderName/$bin/${RUNLIST[$RUN]}"
      echo "alien_cp alien://$PREFIX$bin/$RUN$SUFFIX $bin/${RUNLIST[$RUN]}"
      alien_cp alien://$PREFIX$bin/$RUN$SUFFIX $outputDir/$folderName/$bin/${RUNLIST[$RUN]}
    done
  done
fi

# If run-by-run output not needed, proceed to merge runs, scale, and merge pt-hard bins
if [ "$doRunByRun" = false ]; then

  # (2) Merge output files for each Pt-hard bin
  if [ "$mergeToPtHardBins" = true ]; then
    for bin in $(seq 1 $PtHardBins);
    do
      hadd $outputDir/$folderName/$bin/AnalysisResultsPtHard${bin}.root $outputDir/$folderName/$bin/*/AnalysisResults.root
    done
  fi
  #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  # This is for roughly half of your statistic for unfolding tests
  if [ "$mergeToPtHardBinsSubSample" = true ]; then
    for bin in $(seq 1 $PtHardBins);
    do
      mkdir -p "$outputDir/$folderName/SubSample1/$bin"
      mkdir -p "$outputDir/$folderName/SubSample2/$bin"
      echo "made folders $outputDir/$folderName/SubSample1(2)/$bin/"
      for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
      do
        #Sort the runs into two independent subsets:
        if [ $RUN -lt $NoOfRunsHalf ]; then
          cp "$outputDir/$folderName/$bin/${RUNLIST[$RUN]}/AnalysisResults.root" "$outputDir/$folderName/SubSample1/$bin/Run${RUNLIST[$RUN]}.root"
        else
          cp "$outputDir/$folderName/$bin/${RUNLIST[$RUN]}/AnalysisResults.root" "$outputDir/$folderName/SubSample2/$bin/Run${RUNLIST[$RUN]}.root"
        fi
      done
      hadd $outputDir/$folderName/SubSample1/$bin/AnalysisResultsPtHard${bin}.root $outputDir/$folderName/SubSample1/$bin/Run*.root
      hadd $outputDir/$folderName/SubSample2/$bin/AnalysisResultsPtHard${bin}.root $outputDir/$folderName/SubSample2/$bin/Run*.root
    done
  fi

  # (3) Call python script to re-weight histograms in each Pt-hard bin file
  if [ "$scaleHistograms" = true ]; then
    cd "${outputDir}/${folderName}"
    if [ "$useReferenceFile" = true ]; then
      python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
    else
      python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py
    fi
  fi
  #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #for half of the statistic
  if [ "$scaleHistogramsSubSample" = true -a "$useReferenceFile" = true ]; then
    echo "Scale output with split statistic"
    cd "${outputDir}/${folderName}/SubSample1"
    python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
    cd "${outputDir}/${folderName}/SubSample2"
    python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
  fi

  # (4) Sums all weighted Pt-hard bins into a single final output file.
  if [ "$mergePtHardBins" = true ]; then
    cd "${outputDir}/${folderName}/"
    hadd AnalysisResultsFinal.root */AnalysisResultsPtHard*.root
  fi
  #. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  #for half of the statistic
  if [ "$mergePtHardBinsSubSample" = true ]; then
    cd "${outputDir}/${folderName}/SubSample1/"
    hadd ../AnalysisResultsFinal_SubSample1.root */AnalysisResultsPtHard*.root
    cd "${outputDir}/${folderName}/SubSample2/"
    hadd ../AnalysisResultsFinal_SubSample2.root */AnalysisResultsPtHard*.root
  fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FOR RUN-BY-RUN OUTPUT
else

  # Create dir structure for each run, and merge pt-hard bins for each run
  for bin in $(seq 1 $PtHardBins);
  do
    for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
    do
      if [ ! -d "${outputDir}/${folderName}/${RUNLIST[$RUN]}/${bin}" ]; then
        echo "Creating Folder ${outputDir}/${folderName}/$bin/${RUNLIST[$RUN]}"
        mkdir -p "${outputDir}/${folderName}/${RUNLIST[$RUN]}/$bin"
      fi
      if [ -d "${outputDir}/${folderName}/${RUNLIST[$RUN]}/${bin}" ]; then
        FILE="${outputDir}/${folderName}/${RUNLIST[$RUN]}/${bin}/AnalysisResultsPtHard${bin}.root"
        if [ -f $FILE ]; then
          echo "File $FILE exists."
        else
          echo "File $FILE does not exist."
          echo "copy files into new folder structure"
          cp $outputDir/$folderName/$bin/$RUN/AnalysisResults.root $outputDir/$folderName/$RUN/$bin/AnalysisResultsPtHard$bin.root
        fi
      fi
    done
#if [ -d "$bin" ]; then
#     rm -r $bin
#    fi
  done

  # (3) Call python script to re-weight histograms in each Pt-hard bin file
  if [ "$scaleHistograms" = true ]; then
    for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
    do
      echo "----------------------------------------------------------------------------------------------------"
      echo ""
      echo "Scaling Run ${RUNLIST[$RUN]}"
      echo ""
      echo "----------------------------------------------------------------------------------------------------"
      cd "${outputDir}/${folderName}/${RUNLIST[$RUN]}"
      if [ "$useReferenceFile" = true ]; then
        python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py -f $referenceFile
      else
        python $ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/JetQA/scalePtHardHistos.py
      fi
    done
  fi

  # (4) Sums all weighted Pt-hard bins into a single final output file, for each run.
  if [ "$mergePtHardBins" = true ]; then
    for (( RUN=0; RUN<${NoOfRuns}; RUN++ ));
    do
      cd "${outputDir}/${folderName}/${RUNLIST[$RUN]}"
      echo "----------------------------------------------------------------------------------------------------"
      echo ""
      echo "Merging Run ${RUNLIST[$RUN]}"
      echo ""
      echo "----------------------------------------------------------------------------------------------------"
      hadd AnalysisResultsFinal.root */AnalysisResultsPtHard*.root
    done
  fi

fi


