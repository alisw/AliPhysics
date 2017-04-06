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
#                 /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2012/LHC12i/000193011/lists
###############################################################################
# how to run:  
#
if [ 1 -eq 0 ]; then
    workDir=/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2016/LHC16t
    testDir=/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2016/LHC16t/test2
    codeDir=/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2016/LHC16t/code
    maincodeDir=/lustre/nyx/alice/users/marsland/alice/ali-master/AliPhysics/PWGPP/rawmerge
    cp $maincodeDir/{*.sh,*.C,*.config,*.jdl,*RE*} $codeDir/; rm $codeDir/*~
    cp /u/marsland/PHD/macros/marsland_EbyeRatios/makeOfflineTriggerList.* $codeDir/
    #cd $testDir
    source $codeDir/makeOfflineTriggerList.sh; source $codeDir/makeOfflineTriggerList.config
    
    makeOfflineTriggerListPeriod $workDir/FilteredTrees_lustre_LHC16t.list   
    makeOfflineTriggerListRun  FilteredEvents_267164.list                 ## run a local test in the very beginning
    makeOfflineTriggerListPeriod run267164.list   ## either local or on batchfarm
    copyListsToAlien 2016 LHC16t                              ## run in screen session
    createSubmitJDLlist 2016 LHC16t 50                        ## first make sure that the jdl submit command is ok
    copyRawChunksToLustre 2016 LHC16t /lustre/nyx/alice/alien/triggeredRaw/                        ## run in screen session
fi
###########################################################################################################


if [ "$1" == "-h" ]; then
    echo Usage:     
    echo '(source $ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.sh;  makeOfflineTriggerListPeriod $filteredList)'
    exit;
fi

# below is generic
###########################################################################################################
makeOfflineTriggerListPeriod()
{
  # What it does:
  # Input:
  # 
  # Output: 
  #
  filteredTreeList=$1  
  outbase=`pwd`
  # compile the macro for once
  echo 'gROOT->LoadMacro("$codeDir/makeOfflineTriggerList.C+"); gSystem->Exit(0);' | aliroot -b -l
  # loop over files and run makeOfflineTriggerListRun for each run "FilterEvents_Trees.root" file
  tmprunNumber=0
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
    
    # avoid double counting runnumber and for each run create one list of filtered tree list
    if [ "$runNumber" != "$tmprunNumber" ]; then
      tmprunNumber=$runNumber
    else
      continue
    fi
    runbasedTreeListName=FilteredTrees.list
    $(cat $filteredTreeList | grep $runNumber > $runbasedTreeListName) 
    
    runbasedRawListName=RawFiles.list
    $(alien_find /alice/data/$year/$period/$runNumber/raw .root | grep root | awk '{$0="alien://"$0}1'  > $runbasedRawListName)
    
    runDir=`pwd`/$runNumber/lists;         mkdir -p $runDir; 
    mv `pwd`/$runbasedTreeListName $runDir;  
    mv `pwd`/$runbasedRawListName $runDir;    
    cd $runDir
    # echo "file to be filtered == $runDir/$runbasedTreeListName" 
    # run event merging either locally on farm
    if [ $useBatchFrarm -eq 1 ]; then
      echo "sent job for --> $year $period $runNumber $runDir/$runbasedTreeListName" 
      eval $batchCommand $codeDir/makeOfflineTriggerList.sh  makeOfflineTriggerListRun $runDir/$runbasedTreeListName 
    else
      echo "run locally for --> $year $period $runNumber $runDir/$runbasedTreeListName"
      $codeDir/makeOfflineTriggerList.sh  makeOfflineTriggerListRun $runDir/$runbasedTreeListName  
    fi
    cd $outbase  
  done
}
###########################################################################################################
makeOfflineTriggerListRun()
{
  file=$1
  source $codeDir/makeOfflineTriggerList.config
  echo "*********************************** Input Information *****************************"  
  echo "inputFile      = $file"  
  echo "PWD            = $PWD"
  echo "codeDir        = $codeDir"
  echo "batchCommand   = $batchCommand"
  echo "runType        = $isCosmic"
  echo "useBatchFrarm  = $useBatchFrarm"
  echo "copyTimeout    = $copyTimeoutAlien"
  echo "*********************************** process to be run *****************************" 
  echo aliroot -b -q "${codeDir}/makeOfflineTriggerList.C+(\"${file}\")"
  echo "***********************************************************************************" 
  echo "------------------------- raw event list is being produced ------------------------" 
  # run the macro and filter the output of TTree::Scan 
  aliroot -b -q "${codeDir}/makeOfflineTriggerList.C+(\"${file}\")" &> makeOfflineTriggerList.log 
}
###########################################################################################################
copyListsToAlien()
{
  year=$1
  period=$2  
  alienPrefix="alien:/"
  for file in $(ls $localPrefixRawFilteringDir/$year/$period/000*/lists/RawFiles.list); do
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
    cd $fileDirLocal
    for ifile in $(ls *.*); do 
      echo "file of $run ==  $ifile"
      fileCheck=$(alien_ls $fileDirAlien/$ifile)       
      if [[ "$fileCheck" != *"no such file or directory"* ]]; then 
      echo " ------- $fileDirAlien/$ifile is already copied ------- "
      continue; 
      fi
      alien_cp  $fileDirLocal/$ifile  $alienPrefix/$fileDirAlien
    done 
  done

}
###########################################################################################################
createSubmitJDLlist()
{

  #
  # creates a file containing submit commands for jdl submit 
  # E.g.: 
  #     alien_submit /alice/cern.ch/user/p/pwg_pp/rawmerge/vAN-20170121/mergeOfflineTriggerList.jdl  193094  LHC12i 2016 > submit193094.log &
  #
 
  year=$1
  period=$2
  nfilesPerJob=$3
  alienPrefix="alien:/"

  for run in `alien_ls /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/$year/$period/`; do 
   if [[ "$run" != *"000"* ]]; then continue; fi
   run=${run:3}
   echo  "alien_submit alien:///alice/cern.ch/user/p/pwg_pp/rawmerge/vAN-20170121/mergeOfflineTriggerList.jdl $run $period $year $nfilesPerJob; 2>&1 | tee submit$run.log &"
  done >  submit_$year_$period.sh
  
  # copy the txt file whic contains submitted job list (together with the jdl submit command line)
  fileCheck=$(alien_ls $fileDirAlien/submit_$year_$period.sh)
  if [[ "$fileCheck" != *"no such file or directory"* ]]; then  
    echo " ------- $fileDirAlien/$ifile is already copied renew it ------- "
    alien_rm $fileDirAlien/submit_$year_$period.sh
  fi
  fileDirAlien=$alienPrefixRawFilteringDir/$year/$period
  alien_cp `pwd`/submit_$year_$period.sh $alienPrefix/$fileDirAlien

}
###########################################################################################################
copyRawChunksToLustre()
{
  year=$1
  period=$2 
  outputDir=$3
    
  $(alien_find /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/$year/$period .root | egrep root | egrep rawSelected > tmp.list)
  $(awk '{$0="alien://"$0}1' tmp.list > MergedRawEvents_OnAlien.list); 
  rm tmp.list
  alisync @MergedRawEvents_OnAlien.list -o $outputDir
 
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
###########################################################################################################
ProcessfilterLog(){
    # input parameter  - path to the log prefix
    # log files to be parsed  $logPath/*/lists/makeOfflineTriggerList.log
    # example usage ( source $ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.sh; ProcessfilterLog /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/PWGPP-126/filtering/2016/LHC16t/ )
    # 
    logPath=$1;  
    source $ALICE_PHYSICS/PWGPP/scripts/utilities.sh;     
    egrep KeyValue $logPath/*/lists/makeOfflineTriggerList.log |sed s_":I-KeyValue."_"\t"_ | gawk '{print $1"\t"$2"\t"$3}' >  tmpsummary.log
    echo year/d:period/C:run/d:name/C:key/C:value/I >log.tree
    cat tmpsummary.log | \
	while read line; do
	    guessRunData $line
	    printf "$year\t$period\t$runNumber\t$line\n"
    done >> log.tree
    echo ".L $ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.C" >command.sh
    echo "SummarizeLogs()" >>command.sh
    aliroot -b <command.sh  | tee ProcessfilterLog.log
}


gitRawListFromAlien(){
    #  Example usage 
    #  ( source $ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.sh; gitRawListFromAlien /alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/2015/LHC15n/000244540/ |tee  rawSummary.log )
    # Input:
    #    alien path for run filtered data as stored on 
    # Output:
    #   gidraw.tree               - list with files+gids 
    #   *.download.log            - download log files
    #   info written to std out    - in example above redirected to the rawSummary.log
    #   
    alienPath=$1   
    source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh
    
    # 1.) cache triggered lists locall
    alilog "gitRawList. Step0. DownloadList.Begin"
    for a in `alien_ls $alienPath/lists/filtered*.list`; do echo $alienPath/lists/$a; done  > trigger.list
    alilog "gitRawList. Step0. DownloadList. cat trigger.list" |tee  trigger.download.log
    cat trigger.list
    $ALICE_PHYSICS/PWGPP/QA/scripts/alienSync.sh alienFindCommand="cat trigger.list"
    alilog "gitRawList. Step0. DownloadList.End"
    #
    # 2.) get raw gids - use ls 0stead of find - (it is faster at minimum find is now working very slow)
    #
    alilog "gitRawList. Step1. DownloadLogs.Begin"
    for a in `alien_ls $alienPath`; do echo $alienPath/$a/filter.log; done  > filterlog.list
    alilog "gitRawList. Step0. DownloadList. cat filterlog.list"
    cat filterlog.list
    #$ALICE_ROOT/../src/STEER/Utilities/alienSync.sh alienFindCommand="cat filterlog.list" | tee filterlog.download.log
    $ALICE_PHYSICS/PWGPP/QA/scripts/alienSync.sh alienFindCommand="cat filterlog.list" copyMethod=1  | tee filterlog.download.log
    alilog "gitRawList. Step1. DownloadLogs.End"
    #
    # 2.) Extract event information
    #
    alilog "gitRawList. Step2. ParseLog.Begin"
    echo fname/C:gid/L > gidraw.tree
    find -iname filter.log | xargs cat | grep  AliOfflineTrigger::ExtractSelected: | grep -v N_{trig} | grep -v "Xrd: "| gawk {'printf ("%5s\t%s\n", $3, $5)'}  | sed s_.*raw/__>> gidraw.tree
    echo fname/C:nevents/L > nevents.tree
    find -iname filter.log | xargs cat | grep -v "Xrd: "|  grep  "{trig}" |sed s_.root.*{trig}=_"\ "_g | sed s_.*raw/__ >>nevents.tree
    alilog "gitRawList. Step2. ParseLog.End"
    # 3. Summarize content
    # Number of skipped files
    nFilesSkipped=`find -iname filter.log | xargs cat | grep -c "E-TAlienFile::Open: No more images to try - giving up"`
    nFilesProcessed=`find -iname filter.log | xargs cat |  grep  "I- AliOfflineTrigger::ExtractSelected" | grep -c N_{trig}`
    nFilesEmpty=`find -iname filter.log | xargs cat |  grep  "I- AliOfflineTrigger::ExtractSelected" | grep -c N_{trig}=0`
    nEvents=`find -iname filter.log | xargs cat | grep -v "Xrd: "|  grep  "I- AliOfflineTrigger::ExtractSelected" | grep -v N_{trig}|grep -c root`
    alilog   "gitRawList.nFilesSkipped\t$nFilesSkipped"
    alilog   "gitRawList.nFilesProcessed\t$nFilesProcessed"
    alilog   "gitRawList.nFilesEmpty\t$nFilesEmpty"
    alilog   "gitRawList.nEvents\t$nEvents"
    #
    for a in `find alice -iname "filtered*list" | grep "/lists/" `; do ln -sf $a .  ; done
    ls filtered*list > ref.list
    aliroot -b -q "$ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.C+(\"GetRawSummary\")" | grep KeyValue> GetRawSummary.log

}
  
gitRawListFromAlienPeriod(){
    # this routine for local test purposes in future can be done as part of alien filtering jobs -in merging part
    # period=/alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/2015/LHC15n/
    period=$1
    wdir=`pwd`
    for run in `alien_ls $period | grep 000`; do
	mkdir $wdir/$run;
	alilog "gitRawListFromAlienPeriod.Processing run $run"
        cd $wdir/$run
	( source $ALICE_PHYSICS/PWGPP/rawmerge/makeOfflineTriggerList.sh; gitRawListFromAlien $period/$run/ |tee  rawSummary.log )
    done; 
}



###############################################################################
main()
{
  eval "$@"
}
main "$@"
