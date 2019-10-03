#!/bin/bash

#
# Aim: (Filtered tree --> raw event list) 
# File structures:
#     input: root file  --> root file with "highPt, Laser, V0s, Cosmics" trees:
#                           tree -> Scan() <esdFilePath> <eventNumberinFile> <triggerName> <GID> <timeStamp> <blabla>
#     output:ASCII file -->  <rawFilePath> <eventNumberinFile> <triggerName> <GID> <timeStamp> <blabla>
#     E.g.: 
#     alien:///alice/data/2012/LHC12h/000189122/raw/12000189122032.14.root  1391  CosmicPairs  88801749113   1348661876
#
# Follow the naming convetion:
#            <prefix>/$year/$period/$run/lists
#            E.g: /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/2012/LHC12i/000193011/lists/
#                 /hera/alice/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2012/LHC12i/000193011/lists
###############################################################################
# how to run:  
if [ 1 -eq 0 ]; then
    cd /hera/alice/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2015/LHC15c
    cp /u/marsland/gitTPC_test/alice-tpc-notes/JIRA/PWGPP-7/code/{*.sh,*.C,*.config,*.jdl,*RE*} ./code; rm ./code/*~
    source ./code/makeEventList.sh; source ./code/makeEventList.config
    
    makeEventListRun  FilterEvents_Trees.root                 ## run a local test in the very beginning
    makeEventListPeriod $workDir/FilteredTrees_OnAlien.list   ## either local or on batchfarm
    copyListsToAlien 2012 LHC12i                              ## run in screen session
    createSubmitJDLlist 2012 LHC12i                           ## first make sure that the jdl submit command is ok
    copyRawChunksToHera 2012 LHC12i                           ## run in screen session
fi
###############################################################################

# below is generic
###############################################################################
makeEventListPeriod()
{
  filteredTreeList=$1  
  outbase=`pwd`
  # compile the macro for once
  echo 'gROOT->LoadMacro("$codeDir/makeEventList.C+"); gSystem->Exit(0);' | aliroot -b -l
  # loop over files and run makeEventListRun for each run "FilterEvents_Trees.root" file
  for dataFile in $(cat $filteredTreeList); do
    IFS='/' read -a arr <<< "$dataFile"          # tokenize and create an array another way is arr=(${listName//\// })
    # guess run number etc
    arrLength=$(( ${#arr[*]}-1 ))
    for ((x=0;x<=${arrLength};x++)); do
      field=${arr[${x}]}
      [[ ${field} =~ ^LHC[0-9][0-9][a-z].*$ ]] && period=${field%_*}
      [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && runNumber=${field}
      [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && shortRunNumber=${field#000}
      [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field}
    done
    # run event merging either locally on farm
    runDir=`pwd`/$runNumber/lists; mkdir -p $runDir; cd $runDir
    if [ $useBatchFrarm -eq 1 ]; then
      echo "sent job for --> $year $period $runNumber" 
      eval $batchCommand $codeDir/makeEventList.sh  makeEventListRun $dataFile 
    else
      echo "run locally for --> $year $period $runNumber"
      $codeDir/makeEventList.sh  makeEventListRun $dataFile 
    fi
    cd $outbase  
  done
}
###########################################################################################################
makeEventListRun()
{
  file=$1
  source $codeDir/makeEventList.config
  echo "*********************************** Input Information *****************************"  
  echo "ptMinHighPt    = $ptMinHighPt"
  echo "ptMinV0s       = $ptMinV0s"
  echo "inputFile      = $file"  
  echo "PWD            = $PWD"
  echo "runType        = $isCosmic"
  echo "commonPrefix   = $commonPrefix"
  echo "rawDataPrefix  = $rawDataPrefix" 
  echo "batchCommand   = $batchCommand"
  echo "useBatchFrarm  = $useBatchFrarm"
  echo "alienDirectory = $alienPrefixRawFilteringDir"
  echo "localDirectory = $localPrefixRawFilteringDir"
  echo "codeDir        = $codeDir"
  echo "triggeredRawDir= $heraTriggeredRawDir"
  echo "copyTimeout    = $copyTimeoutAlien"
  echo "*********************************** process to be run *****************************" 
  echo aliroot -b -q "${codeDir}/makeEventList.C+(\"${file}\",${ptMinHighPt},${ptMinV0s})"
  echo "***********************************************************************************" 
  echo "------------------------- raw event list is being produced ------------------------" 
  # run the macro and filter the output of TTree::Scan 
  aliroot -b -q "${codeDir}/makeEventList.C+(\"${file}\",${ptMinHighPt},${ptMinV0s})" &> makeEventList.log 
}
###########################################################################################################
copyListsToAlien()
{
  year=$1
  period=$2  
  alienPrefix="alien:/"
  for file in $(ls $localPrefixRawFilteringDir/$year/$period/*/lists/event.list); do
    # estimate run number from path
    IFS='/' read -a arr <<< "$file" # tokenize and create an array another way is arr=(${listName//\// })
    arrLength=${#arr[@]}; 
    for ((x=0;x<=${arrLength};x++)); do
      field=${arr[${x}]}
      [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && run=${field}
    done
    fileDirAlien=$alienPrefixRawFilteringDir/$year/$period/$run/lists
    fileDirLocal=$localPrefixRawFilteringDir/$year/$period/$run/lists
    echo "fileDirAlien   = " $fileDirAlien
    echo "fileDirLocal   = " $fileDirLocal
    alien_mkdir -p $fileDirAlien
    for ilist in $(ls $fileDirLocal/*.list); do 
      echo $ilist
      alien_cp  $ilist  $alienPrefix/$fileDirAlien
    done 
  done

}
###########################################################################################################
createSubmitJDLlist(){

  #
  # creates a file containing submit commands for jdl submit 
  # E.g.: 
  #     alien_submit /alice/cern.ch/user/p/pwg_pp/rawmerge/v20150512/rawmerge.jdl  193094  LHC12i 2012 > submit193094.log &
  #
 
  year=$1
  period=$2
  alienPrefix="alien:/"

  for run in `alien_ls /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/$year/$period/  | grep 000 | sed s_000__ `; do 
   echo  "alien_submit alien:///alice/cern.ch/user/p/pwg_pp/rawmerge/v20150512/rawmerge.jdl $run  $period $year; 2>&1 | tee submit$run.log &"
  done >  submit_$year_$period.sh
  
  # copy the txt file whic contains submitted job list (together with the jdl submit command line)
  fileDirAlien=$alienPrefixRawFilteringDir/$year/$period
  alien_cp `pwd`/submit_$year_$period.sh $alienPrefix/$fileDirAlien

}
###########################################################################################################
copyRawChunksToHera()
{
  year=$1
  period=$2    
  findCommand="alien_find /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/$year/$period root"
  $ALICE_PHYSICS/../src/PWGPP/QA/scripts/alienSync.sh alienFindCommand=$findCommand localPathPrefix=$heraTriggeredRawDir copyTimeout=$copyTimeoutAlien
}
###########################################################################################################
DirCheckCreate()
{
  dirName=$1
  if [ -d "$dirName" ]; then
  echo " !!! Attention !!! " $dirName "  exists already. Delete and recreate  "
  rm -rf $dirName 
  fi
  mkdir $dirName
}
###############################################################################
main()
{
  eval "$@"
}
main "$@"
