#!/bin/bash
# this script runs the CPass0/CPass1 train
# produced OCDB updates are local

main()
{
  #run in proper mode depending on the selection
  runMode=$1
  if [[ $# -gt 0 ]]; then
    echo "# $0 $*"
  fi
  umask 0002
  shift
  case $runMode in
    "CPass0") goCPass0 "$@";;
    "CPass1") goCPass1 "$@";;
    "ConfOCDB") goConfigureOCDBforCPass1 "$@";;
    "MergeCPass0") goMergeCPass0 "$@";;
    "MergeCPass1") goMergeCPass1 "$@";;
    "CreateQAplots") goCreateQAplots "$@";;
    "WaitForOutput") goWaitForOutput "$@";;
    "makeSummary") goMakeSummary "$@";;
    "submit") goSubmit "$@";;
    "submitQA") goSubmitQA "$@";;
    "test") goTest "$@";;
    "generateMakeflow") goGenerateMakeflow "$@";;
    "makeflow") goMakeflow "$@";;
    #and the default: unless sourced print some info
    *) if [[ ! "$0" =~ "bash" ]]; then
         echo " batch:    $0 \"submit\" inputList productionID [configFile=benchmark.config] [runNumber]"
         echo " makeflow: $0 \"makeflow\" inputList productionID [configFile=benchmark.config] [runNumber]"
       fi
    ;;
  esac
}

goCPass0()
{
  umask 0002
  
  configFile=$5
  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} && -z ${alirootEnv} ]] && source $setupAliROOTenvInCurrentShell

  targetDirectory=$1
  inputList=$2
  nEvents=$3
  ocdbPath=$4
  configFile=$5
  runNumber=$6
  jobindex=$7

  #use the jobindex only if set and non-negative
  if [[ -z $jobindex || $jobindex -lt 0 ]]; then
    [[ -n "$LSB_JOBINDEX" ]] && jobindex=$LSB_JOBINDEX
    [[ -n "$SGE_TASK_ID" ]] && jobindex=$SGE_TASK_ID
  fi

  [[ ! -f ${inputList} && -z ${pretend} ]] && echo "input file $inputList not found, exiting..." && exit 1
  if [[ "${inputList}" =~ \.root$ ]]; then
    infile=$inputList
  else
    infile=`sed -ne "${jobindex}p" $inputList`
  fi
  
  chunkName=${infile##*/}
  outputDir=${targetDirectory}/${jobindex}
  mkdir -p $outputDir
  [[ ! -d $outputDir ]] && echo "cannot make $outputDir" && exit 1
  
  runpath=${PWD}/rundir_cpass0_${runNumber}_${jobindex}
  [[ -z $commonOutputPath ]] && commonOutputPath=$PWD
  [[ $reconstructInTemporaryDir -eq 1 && -n $TMPDIR ]] && runpath=$TMPDIR
  [[ $reconstructInTemporaryDir -eq 1 && -z $TMPDIR ]] && runpath=$(mktemp -d)

  mkdir -p $runpath
  [[ ! -d ${runpath} ]] && echo "cannot make runpath ${runpath}" && exit 1
  cd $runpath

  logOutputDir=$runpath
  [[ -n $logToFinalDestination ]] && logOutputDir=$outputDir
  [[ -z $dontRedirectStdOutToLog ]] && exec 1> $logOutputDir/stdout
  [[ -z $dontRedirectStdOutToLog ]] && exec 2> $logOutputDir/stderr
  echo "$0 $*"

  calibDoneFile="$commonOutputPath/cpass0.job${jobindex}.run${runNumber}.done"

  echo "#####################"
  echo CPass0:
  echo JOB setup
  echo nEvents            $nEvents
  echo runNumber          $runNumber
  echo ocdbPath           $ocdbPath
  echo infile             $infile
  echo chunkName          $chunkName
  echo jobindex           $jobindex
  echo recoTriggerOptions $recoTriggerOptions
  echo targetDirectory    $targetDirectory
  echo commonOutputPath         $commonOutputPath
  echo calibDoneFile      $calibDoneFile
  echo runpath            $runpath  
  echo outputDir          $outputDir
  echo ALICE_ROOT         $ALICE_ROOT
  echo PWD                $PWD
  echo "########## ###########"

  alirootInfo > ALICE_ROOT_svn.log

  filesCPass0=( 
               "$commonOutputPath/runCPass0.sh"
               "$commonOutputPath/recCPass0.C"
               "$commonOutputPath/runCalibTrain.C"
               "$commonOutputPath/localOCDBaccessConfig.C"
               "$commonOutputPath/OCDB.root"
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/runCPass0.sh"
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/recCPass0.C" 
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/runCalibTrain.C"
  )

  for file in ${filesCPass0[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done

  ln -s $infile $runpath/$chunkName

  echo "this directory ($PWD) contents:"
  ls -lh
  echo
  chmod u+x runCPass0.sh

  if [[ -n $postSetUpActionCPass0 ]]; then
    echo "running $postSetUpActionCPass0"
    eval $postSetUpActionCPass0
  fi

  #run CPass0
  echo "$runpath/runCPass0.sh $infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
  if [[ -n $pretend ]]; then
    touch AliESDfriends_v1.root
    touch rec.log
    touch calib.log
  else
    ./runCPass0.sh "$infile" "$nEvents" "$runNumber" "$ocdbPath" "$recoTriggerOptions"
  fi
  
  #move stuff to final destination
  echo "this directory ($PWD) contents:"
  ls -lh
  echo

  echo rm -f ./$chunkName
  rm -f ./$chunkName
  echo "cp --recursive $runpath/* ${outputDir}"
  cp --recursive $runpath/* $outputDir
  echo
  
  #validate CPass0
  cd ${outputDir}
  touch ${calibDoneFile}
  [[ -f AliESDfriends_v1.root ]] && echo "calibfile ${outputDir}/AliESDfriends_v1.root" > ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}

  rm -rf ${runpath}
}

goCPass1()
{
  umask 0002
  
  configFile=$5
  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} && -z ${alirootEnv} ]] && source $setupAliROOTenvInCurrentShell

  targetDirectory=$1
  inputList=$2
  nEvents=$3
  ocdbPath=$4
  configFile=$5
  runNumber=$6
  jobindex=$7

  #use the jobindex only if set and non-negative
  if [[ -z $jobindex || $jobindex -lt 0 ]]; then
    [[ -n "$LSB_JOBINDEX" ]] && jobindex=$LSB_JOBINDEX
    [[ -n "$SGE_TASK_ID" ]] && jobindex=$SGE_TASK_ID
  fi

  [[ ! -f ${inputList} && -z ${pretend} ]] && echo "input file $inputList not found, exiting..." && exit 1
  if [[ "${inputList}" =~ \.root$ ]]; then
    infile=$inputList
  else
    infile=`sed -ne "${jobindex}p" $inputList`
  fi
  
  chunkName=${infile##*/}
  outputDir=${targetDirectory}/${jobindex}
  mkdir -p $outputDir
  [[ ! -d $outputDir ]] && echo "cannot make $outputDir" && exit 1
  
  runpath=${PWD}/rundir_cpass1_${runNumber}_${jobindex}
  [[ -z $commonOutputPath ]] && commonOutputPath=$PWD
  [[ $reconstructInTemporaryDir -eq 1 && -n $TMPDIR ]] && runpath=$TMPDIR
  [[ $reconstructInTemporaryDir -eq 1 && -z $TMPDIR ]] && runpath=$(mktemp -d)

  #init the running path
  mkdir -p $runpath
  [[ ! -d ${runpath} ]] && echo "cannot make runpath ${runpath}" && exit 1
  cd $runpath


  calibDoneFile="${commonOutputPath}/cpass1.job${jobindex}.run${runNumber}.done"

  logOutputDir=$runpath
  [[ -n $logToFinalDestination ]] && logOutputDir=$outputDir
  [[ -z $dontRedirectStdOutToLog ]] && exec 1> $logOutputDir/stdout
  [[ -z $dontRedirectStdOutToLog ]] && exec 2> $logOutputDir/stderr
  echo "$0 $*"

  echo "#####################"
  echo CPass1:
  echo JOB setup
  echo nEvents            $nEvents
  echo runNumber          $runNumber
  echo ocdbPath           $ocdbPath
  echo infile             $infile
  echo chunkName          $chunkName
  echo jobindex           $jobindex
  echo recoTriggerOptions $recoTriggerOptions
  echo targetDirectory    $targetDirectory
  echo commonOutputPath         $commonOutputPath
  echo calibDoneFile      $calibDoneFile
  echo runpath            $runpath  
  echo outputDir          $outputDir
  echo ALICE_ROOT         $ALICE_ROOT
  echo PWD                $PWD
  echo "########## ###########"

  alirootInfo > ALICE_ROOT_svn.log

  filesCPass1=( 
               "$commonOutputPath/runCPass1.sh"
               "$commonOutputPath/recCPass1.C"
               "$commonOutputPath/recCPass1_OuterDet.C"
               "$commonOutputPath/runCalibTrain.C"
               "$commonOutputPath/QAtrain_duo.C"
               "$commonOutputPath/localOCDBaccessConfig.C"
               "$commonOutputPath/cpass0.localOCDB.${runNumber}.tgz"
               "$commonOutputPath/OCDB.root"
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/runCPass1.sh"
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/recCPass1.C" 
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/recCPass1_OuterDet.C" 
               "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/runCalibTrain.C"
               "$ALICE_ROOT/ANALYSIS/macros/QAtrain_duo.C"
  )

  for file in ${filesCPass1[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done

  ln -s $infile $runpath/$chunkName

  echo "this directory ($PWD) contents:"
  ls -lh
  echo

  if [[ -n $postSetUpActionCPass1 ]]; then
    echo "running $postSetUpActionCPass1"
    eval $postSetUpActionCPass1
    echo
  fi

  #configure local OCDB storage from CPass0 (creates the localOCDBaccessConfig.C script)
  if [[ -f cpass0.localOCDB.${runNumber}.tgz ]]; then
   echo goConfigureOCDBforCPass1 "cpass0.localOCDB.${runNumber}.tgz"
   goConfigureOCDBforCPass1 "cpass0.localOCDB.${runNumber}.tgz"
 else
   echo "WARNING: file cpass0.localOCDB.${runNumber}.tgz not found!"
 fi

  #run CPass1
  chmod u+x runCPass1.sh
  echo "$runpath/runCPass1.sh $infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
  if [[ -n $pretend ]]; then
    touch AliESDfriends_v1.root
    touch QAresults_Barrel.root
    touch QAresults_Outer.root
    touch rec.log
    touch calib.log
    touch qa.log
  else
    ./runCPass1.sh "$infile" "$nEvents" "$runNumber" "$ocdbPath" "$recoTriggerOptions"
  fi
  
  #move stuff to final destination
  echo "this directory ($PWD) contents:"
  ls
  echo

  echo rm -f ./$chunkName
  rm -f ./$chunkName
  echo "cp --recursive ${runpath}/* ${outputDir}"
  cp --recursive ${runpath}/* ${outputDir}
  echo

  #validate CPass1
  cd ${outputDir}
  touch ${calibDoneFile}
  [[ -f AliESDfriends_v1.root ]] && echo "calibfile ${outputDir}/AliESDfriends_v1.root" > ${calibDoneFile}
  [[ -f QAresults_Barrel.root ]] && echo "qafile ${outputDir}/QAresults_Barrel.root" >> ${calibDoneFile}
  [[ -f QAresults_Outer.root ]] && echo "qafile ${outputDir}/QAresults_Outer.root" >> ${calibDoneFile}
  echo "dir ${outputDir}" >> ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}
  
  rm -rf ${runpath}
}


goMergeCPass0()
{
  umask 0002
  #
  # find the output files and merge them
  #

  outputDir=$1
  defaultOCDB=$2
  configFile=$3
  runNumber=$4
  calibrationFilesToMergeExternal=$5

  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} && -z ${alirootEnv} ]] && source $setupAliROOTenvInCurrentShell

  runpath=${PWD}/rundir_cpass0_Merge_${runNumber}
  [[ -z $commonOutputPath ]] && commonOutputPath=$PWD
  [[ $reconstructInTemporaryDir -eq 1 && -n $TMPDIR ]] && runpath=$TMPDIR
  [[ $reconstructInTemporaryDir -eq 1 && -z $TMPDIR ]] && runpath=$(mktemp -d)

  mkdir -p $runpath
  [[ ! -d $runpath ]] && echo "not able to make the runpath $runpath" && exit 1
  cd $runpath

  logOutputDir=$runpath
  [[ -n $logToFinalDestination ]] && logOutputDir=$outputDir
  [[ -z $dontRedirectStdOutToLog ]] && exec 2>&1 > $logOutputDir/mergeMakeOCDB.log
  echo "$0 $*"

  calibrationFilesToMerge=$calibrationFilesToMergeExternal
  [[ -z $calibrationFilesToMerge ]] && calibrationFilesToMerge="calibrationFilesToMerge.list"
  calibrationOutputFileName="AliESDfriends_v1.root"
  mergingScript="mergeMakeOCDB.byComponent.sh"
  qaFilesToMerge="qaFilesToMerge.list"
  qaOutputFileName="QAresults*.root"
  qaMergedOutputFileName="QAresults_merged.root"

  echo goMergeCPass0 SETUP:
  echo runNumber=$runNumber
  echo outputDir=$outputDir
  echo defaultOCDB=$defaultOCDB
  echo calibrationFilesToMerge=$calibrationFilesToMerge
  echo calibrationOutputFileName=$calibrationOutputFileName
  echo mergingScript=$mergingScript
  
  # copy files in case they are not already there
  filesMergeCPass0=(
                    "$commonOutputPath/${calibrationFilesToMerge}"
                    "$commonOutputPath/OCDB.root"
                    "$commonOutputPath/localOCDBaccessConfig.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeMakeOCDB.byComponent.sh"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeByComponent.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/merge.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/mergeMakeOCDB.sh"
  )
  for file in ${filesMergeCPass0[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done
  
  alirootInfo > ALICE_ROOT_svn.log

  #
  ls -lh

  #merge calibration
  chmod u+x $mergingScript  
  mkdir -p ./OCDB
  if [[ -z ${calibrationFilesToMergeExternal} ]]; then
    echo "find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge"
    find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge
  fi
  
  echo "$mergingScript $calibrationFilesToMerge ${runNumber} local://./OCDB $defaultOCDB"
  if [[ -n $pretend ]]; then
    touch CalibObjects.root
    touch ocdb.log
    touch merge.log
    mkdir -p ./OCDB/someDetector/
    mkdir -p ./OCDB/otherDetector/
    touch ./OCDB/someDetector/someCalibObject_0-999999_cpass0.root
    touch ./OCDB/otherDetector/otherCalibObject_0-999999_cpass0.root
  else
    ./$mergingScript $calibrationFilesToMerge ${runNumber} "local://./OCDB" $defaultOCDB
  fi

  ### produce the output
  #tar the produced OCDB for reuse
  tar czf $commonOutputPath/cpass0.localOCDB.${runNumber}.tgz ./OCDB

  ls -ltrh

  #copy all to output dir
  cp --recursive ${runpath}/* $outputDir
  
  #validate merging cpass0
  cd ${outputDir}
  calibDoneFile="${commonOutputPath}/merge.cpass0.run${runNumber}.done"
  touch ${calibDoneFile}
  [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" > ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}

  rm -rf ${runpath}
}

goMergeCPass1()
{
  umask 0002
  #
  # find the output files and merge them
  #

  outputDir=$1
  defaultOCDB=$2
  configFile=$3
  runNumber=$4
  calibrationFilesToMergeExternal=$5
  qaFilesToMergeExternal=$6

  #clean up first:
  rm -f $outputDir/*.log
  rm -f $outputDir/*.root
  rm -f $outputDir/*done

  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} && -z ${alirootEnv} ]] && source $setupAliROOTenvInCurrentShell

  runpath=${PWD}/rundir_cpass1_Merge_${runNumber}
  [[ -z $commonOutputPath ]] && commonOutputPath=$PWD
  [[ $reconstructInTemporaryDir -eq 1 && -n $TMPDIR ]] && runpath=$TMPDIR
  [[ $reconstructInTemporaryDir -eq 1 && -z $TMPDIR ]] && runpath=$(mktemp -d)

  mkdir -p $runpath
  [[ ! -d $runpath ]] && echo "not able to make the runpath $runpath" && exit 1
  cd $runpath

  logOutputDir=$runpath
  [[ -n $logToFinalDestination ]] && logOutputDir=$outputDir
  [[ -z $dontRedirectStdOutToLog ]] && exec 2>&1 > $logOutputDir/mergeMakeOCDB.log
  echo "$0 $*"

  calibrationFilesToMerge=$calibrationFilesToMergeExternal
  [[ -z $calibrationFilesToMerge ]] && calibrationFilesToMerge="calibrationFilesToMerge.list"
  calibrationOutputFileName="AliESDfriends_v1.root"
  mergingScript="mergeMakeOCDB.byComponent.sh"
  qaFilesToMerge=$qaFilesToMergeExternal
  [[ -z $qaFilesToMerge ]] && qaFilesToMerge="qaFilesToMerge.list"
  qaOutputFileName="QAresults*.root"
  qaMergedOutputFileName="QAresults_merged.root"

  echo goMergeCPass1 SETUP:
  echo runNumber=$runNumber
  echo outputDir=$outputDir
  echo defaultOCDB=$defaultOCDB
  echo calibrationFilesToMerge=filesToMerge.list
  echo calibrationOutputFileName=$calibrationOutputFileName
  echo mergingScript=$mergingScript
  
  # copy files in case they are not already there
  filesMergeCPass1=(
                    "$commonOutputPath/${calibrationFilesToMerge}"
                    "$commonOutputPath/${qaFilesToMerge}"
                    "$commonOutputPath/OCDB.root"
                    "$commonOutputPath/localOCDBaccessConfig.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.byComponent.sh"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeByComponent.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/makeOCDB.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/merge.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.sh"
  )
  for file in ${filesMergeCPass1[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done

  alirootInfo > ALICE_ROOT_svn.log

  #
  ls -lh

  #merge calibration
  chmod u+x $mergingScript  
  mkdir -p OCDB
  if [[ -z ${calibrationFilesToMergeExternal} ]]; then
    echo "find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge"
    find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge
  fi
  
  echo "$mergingScript $calibrationFilesToMerge ${runNumber} local://./OCDB $defaultOCDB"
  if [[ -n $pretend ]]; then
    touch CalibObjects.root
    touch ocdb.log
    touch merge.log
  else
    ./$mergingScript $calibrationFilesToMerge ${runNumber} "local://./OCDB" $defaultOCDB
  fi

  tar czf localCPass1_${runNumber}.tgz ./OCDB

  #merge QA
  [[ -n ${AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF
  [[ -n ${AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF

  if [[ -z $qaFilesToMergeExternal ]]; then
    echo "find $outputDir -name $qaOutputFileName > $qaFilesToMerge"
    find $outputDir -name $qaOutputFileName > $qaFilesToMerge
  fi
  
  echo aliroot -l -b -q "merge.C(\"$qaFilesToMerge\",\"\",kFALSE,\"$qaMergedOutputFileName\")"
  if [[ -n $pretend ]]; then
    touch $qaMergedOutputFileName
    touch merge.log
  else
    aliroot -l -b -q "merge.C(\"$qaFilesToMerge\",\"\",kFALSE,\"$qaMergedOutputFileName\")"
  fi

  ls -ltrh

  #copy all to output dir
  cp --recursive ${runpath}/* ${outputDir}
  
  #validate merge cpass1
  cd ${outputDir}
  calibDoneFile="${commonOutputPath}/merge.cpass1.run${runNumber}.done"
  touch ${calibDoneFile}
  [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" > ${calibDoneFile}
  [[ -f $qaMergedOutputFileName ]] && echo "qafile $outputDir/$qaMergedOutputFileName" >> ${calibDoneFile}
  echo "dir ${outputDir}" >> ${calibDoneFile}
  summarizeLogs >>  ${calibDoneFile}

  rm -rf ${runpath}
}

goMakeSummary()
{
  configFile=$1
  source $configFile
  
  # log filtering, script needs to take the base dir as argument
  if [[ -x $logFilteringScript ]]; then
    commonOutputPath=${baseOutputDirectory}/${productionID}
    ${logFilteringScript} $commonOutputPath
  fi

  awk 'BEGIN {nFiles=0;} /^calibfile/ {nFiles++;} END {print     "cpass0 produced "nFiles" calib files";}' cpass0.job*done
  awk 'BEGIN {nOK=0; nBAD=0;} /rec.*log OK/ {nOK++;} /rec.*log BAD/ {nBAD++;} END {print     "cpass0 reco:  OK: "nOK" BAD: "nBAD;}' cpass0.job*done
  awk 'BEGIN {nOK=0; nBAD=0;} /calib.*log OK/ {nOK++;} /calib.*log BAD/ {nBAD++;} END {print "cpass0 calib: OK: "nOK" BAD: "nBAD;}' cpass0.job*done
  
  awk 'BEGIN {nOK=0; nBAD=0;} /merge.*log OK/ {nOK++;} /merge.*log BAD/ {nBAD++;} END {print "cpass0 merge: OK: "nOK" BAD: "nBAD;}' merge.cpass0*done
  awk 'BEGIN {nOK=0; nBAD=0;} /ocdb.*log OK/ {nOK++;} /ocdb.*log BAD/ {nBAD++;} END {print   "cpass0 OCDB:  OK: "nOK" BAD: "nBAD;}' merge.cpass0*done
  
  awk 'BEGIN {nFiles=0;} /^calibfile/ {nFiles++;} END {print     "cpass1 produced "nFiles" calib files";}' cpass1.job*done
  awk 'BEGIN {nOK=0; nBAD=0;} /rec.*log OK/ {nOK++;} /rec.*log BAD/ {nBAD++;} END {print     "cpass1 reco:  OK: "nOK" BAD: "nBAD;}' cpass1.job*done
  awk 'BEGIN {nOK=0; nBAD=0;} /calib.*log OK/ {nOK++;} /calib.*log BAD/ {nBAD++;} END {print "cpass1 calib: OK: "nOK" BAD: "nBAD;}' cpass1.job*done

  awk 'BEGIN {nOK=0; nBAD=0;} /merge.*log OK/ {nOK++;} /merge.*log BAD/ {nBAD++;} END {print "cpass1 merge: OK: "nOK" BAD: "nBAD;}' merge.cpass1*done
  awk 'BEGIN {nOK=0; nBAD=0;} /ocdb.*log OK/ {nOK++;} /ocdb.*log BAD/ {nBAD++;} END {print   "cpass1 OCDB:  OK: "nOK" BAD: "nBAD;}' merge.cpass1*done
  
  #if set email the summary
  [[ -n $mailSummaryTo ]] && cat $log | mail -s "benchmark $productionID done" $mailSummaryTo

  return 0
}

goMakeflow()
{
  #generate the makeflow file and run
  inputFileList=$1
  productionID=$2
  configFile=$3
  runNumber=$4

  [[ -z ${configFile} ]] && configFile="benchmark.config"
  [[ ! -f ${configFile} ]] && echo "no config file found (${configFile})" && return 1
  source $configFile

  source $configFile
  goGenerateMakeflow "$@" > benchmark.makeflow
  makeflow ${makeflowOptions} benchmark.makeflow
}

goGenerateMakeflow()
{
  #generate the makeflow file
  inputFileList=$1
  productionID=$2
  configFile=$3
  runNumber=$4

  [[ -z ${configFile} ]] && configFile="benchmark.config"
  [[ ! -f ${configFile} ]] && echo "no config file found (${configFile})" && return 1
  source $configFile

  commonOutputPath=${baseOutputDirectory}/${productionID}

  #these files will be made a dependency - will be copied to the working dir of the jobs
  declare -a copyFiles
  inputFiles=(
              "OCDB.root"
              "localOCDBaccessConfig.C"
  )
  for file in ${inputFiles[*]}; do
    [[ -f ${file} ]] && copyFiles+=("${file}")
  done

  #create the makeflow file
  declare -a arr_cpass1_final
  declare -a arr_cpass1_QA_final
  listOfRuns=${runNumber}
  [[ -z ${runNumber} ]] && listOfRuns=($(while read x; do guessRunNumber $x; done < ${inputFileList} | sort | uniq))
  runindex=0
  for runNumber in ${listOfRuns[*]}; do
    [[ -z $runNumber ]] && continue
    [[ ! ${runNumber} =~ ^[0-9]*[0-9]$ ]] && continue
    jobindex=0

    unset arr_cpass0_outputs
    unset arr_cpass1_outputs
    declare -a arr_cpass0_outputs
    declare -a arr_cpass1_outputs
    unset arr_cpass0_outputs

    while read inputFile; do
      currentDefaultOCDB=${defaultOCDB}
      [[ ${autoOCDB} -ne 0 ]] && currentDefaultOCDB=$(setYear $inputFile $defaultOCDB)

      #CPass0
      arr_cpass0_outputs[$jobindex]="cpass0.job${jobindex}.run${runNumber}.done"
      echo "${arr_cpass0_outputs[$jobindex]}: benchmark.sh ${configFile} ${copyFiles[@]}"
      echo " ${alirootEnv} ./benchmark.sh CPass0 ${commonOutputPath}/cpass0/000${runNumber} $inputFile $nEvents $currentDefaultOCDB $configFile $runNumber $jobindex"
      echo

      #CPass1
      arr_cpass1_outputs[$jobindex]="cpass1.job${jobindex}.run${runNumber}.done"
      echo "${arr_cpass1_outputs[$jobindex]} : benchmark.sh ${configFile} cpass0.localOCDB.${runNumber}.tgz ${copyFiles[@]}"
      echo " ${alirootEnv} ./benchmark.sh CPass1 ${commonOutputPath}/cpass1/000${runNumber} $inputFile $nEvents $currentDefaultOCDB $configFile $runNumber $jobindex"
      echo
      ((jobindex++))

    done< <(grep "/000$runNumber/" $inputFileList)

    #CPass0 list of Calib files to merge
    echo "cpass0.calib.run${runNumber}.list: ${arr_cpass0_outputs[*]}"
    echo "  awk '/^calibfile / {print "'\$2'"}' ${arr_cpass0_outputs[*]} > cpass0.calib.run${runNumber}.list"
    echo

    #CPass1 list of Calib/QA files to merge
    echo "cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.list: ${arr_cpass1_outputs[*]}"
    echo "  awk '/^calibfile / {print "'\$'"2}' ${arr_cpass1_outputs[*]} > cpass1.calib.run${runNumber}.list;  awk '/^qafile / {print "'\$'"2}' ${arr_cpass1_outputs[*]} > cpass1.QA.run${runNumber}.list"
    echo

    #CPass0 merging
    arr_cpass0_final[$runindex]="merge.cpass0.run${runNumber}.done"
    echo "cpass0.localOCDB.${runNumber}.tgz ${arr_cpass0_final[$runindex]}: cpass0.calib.run${runNumber}.list benchmark.sh ${configFile} ${copyFiles[@]}"
    echo " ${alirootEnv} ./benchmark.sh MergeCPass0 ${commonOutputPath}/cpass0/000${runNumber} $currentDefaultOCDB ${configFile} $runNumber cpass0.run${runNumber}.list"
    echo

    #CPass1 Calib/QA merging
    arr_cpass1_final[$runindex]="merge.cpass1.run${runNumber}.done"
    echo "${arr_cpass1_final[$runindex]}: cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.list benchmark.sh ${configFile} ${copyFiles[@]}"
    echo " ${alirootEnv} ./benchmark.sh MergeCPass1 ${commonOutputPath}/cpass1/000${runNumber} $currentDefaultOCDB ${configFile} $runNumber cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.list"
    echo
    ((runindex++))
  done

  #CPass1 list of final Calib/QA files
  echo "cpass1.QA.list cpass1.calib.list: ${arr_cpass1_final[*]}"
  echo " awk '/^calibfile / {print "'\$'"2}' ${arr_cpass1_final[*]} > cpass1.calib.list; awk '/^qafile / {print "'\$'"2}' ${arr_cpass1_final[*]} > cpass1.QA.list"
  echo

  #Summary
  echo "summary.log : ${arr_cpass0_outputs[*]} ${arr_cpass1_outputs[*]}  ${arr_cpass1_final[*]}  ${arr_cpass0_final[*]} benchmark.sh ${configFile}"
  echo " LOCAL ./benchmark.sh makeSummary ${configFile} |tee summary.log"
}

goCreateQAplots()
{
  umask 0002
  [[ $# -lt 5 ]] && echo "goCreateQAplots productionID pass outputDir qaFilesDirectory qaPlotScript" && exit 1
  productionID=$1
  pass=$2
  outputDir=$3
  qaFilesDirectory=$4
  qaPlotsScript=$5
  configFile=$6

  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} ]] && source $setupAliROOTenvInCurrentShell

  runpath=${PWD}/rundir_cpass0_Merge_${runNumber}
  [[ -z $commonOutputPath ]] && commonOutputPath=$PWD
  [[ $reconstructInTemporaryDir -eq 1 && -n $TMPDIR ]] && runpath=$TMPDIR
  [[ $reconstructInTemporaryDir -eq 1 && -z $TMPDIR ]] && runpath=$(mktemp -d)

  mkdir -p $runpath
  [[ ! -d $runpath ]] && echo "not able to make the runpath $runpath" && exit 1
  cd $runpath

  [[ -z $logOutputDir ]] && logOutputDir=$runpath
  [[ -z $dontRedirectStdOutToLog ]] && exec 2>&1 > $logOutputDir/makeQAplots.log
  echo "$0 $*"

  [[ -z "$qaPlotsScript" ]] && echo "qaPlotsScript not defined"&&exit 1
  
  mergedQAfileList=$outputDir/mergedQAfiles.list
  echo "MakeListOfQAresults $qaFilesDirectory QAresults.root | grep CPass1 > $mergedQAfileList"
  MakeListOfQAresults $qaFilesDirectory QAresults_merged.root | grep CPass1 |tee $mergedQAfileList
  echo $qaPlotsScript "$productionID" "cpass1" $mergedQAfileList $outputDir
  $qaPlotsScript "$productionID" "cpass1" $mergedQAfileList $outputDir

  mkdir -p $outputDir
  [[ ! -d $outputDir ]] && echo "cannot make the output dir $outputDir" && exit 1
  mv -f $runpath/* $outputDir
  rm -rf $runpath
}

goWaitForOutput()
{
  umask 0002
  [[ $# -lt 3 ]] && echo "goWaitForOutput() wrong number of arguments, exiting.." && exit 1
  echo searchPath=$1
  echo fileName=$2
  echo numberOfFiles=$3
  echo maxSecondsToWait=$4
  searchPath=$1
  fileName=$2
  numberOfFiles=$3
  maxSecondsToWait=$4
  extraFindOptions=$5
  echo "command to be executed: find $searchPath -name "$fileName" ${extraFindOptions}"
  [[ -z "$maxSecondsToWait" ]] && maxSecondsToWait=$(( 3600*12 ))
  while sleep 60; do
    n=$(find $searchPath -name "$fileName" ${extraFindOptions}| wc -l)
    [[ $n -gt 0 ]] && echo "found $n X $fileName"
    [[ $n -ge $numberOfFiles ]] && break
    [[ $SECONDS -gt $maxSecondsToWait ]] && break
  done
  echo "DONE! exiting..."
}

submit()
{
  umask 0002
  [[ $# -ne 6 ]] && echo "6 args needed, you supplied $#" && exit 1
  JobID=$1
  startID=$2
  endID=$3
  waitForJOBID=$4
  command=$5
  commandArgs=$6

  newFarm=$(which qsub|grep "^/usr/bin/qsub")
  
  batchSystem="SGE"

  if [[ -z "$newFarm" ]]
  then
    #old LSF
    # submit it (as job array)
    nFiles=$(( $endID-$startID+1 ))
    while [ $startID -le $nFiles  ] ; do
      if [ `expr $nFiles - $startID` -gt 999 ] ; then 
        endID=`expr $startID + 999`
      else
        endID=$nFiles
      fi      
      if [[ -z "$waitForJOBID" ]]; then
        echo $batchCommand -J "$JobID[$startID-$endID]" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
        $batchCommand -J "$JobID[$startID-$endID]" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
      else
        echo $batchCommand -J "$JobID[$startID-$endID]" -w "ended($waitForJOBID)" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
        $batchCommand -J "$JobID[$startID-$endID]" -w "ended($waitForJOBID)" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
      fi
      startID=`expr $endID + 1`
    done
  else 
    #new SGE farm
    if [[ -z "$waitForJOBID" ]]; then
      echo $batchCommand -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
      $batchCommand -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
    else
      echo $batchCommand -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -hold_jid "$waitForJOBID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
      $batchCommand -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -hold_jid "$waitForJOBID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
    fi
  fi
}

goTest()
{
  umask 0002
  exec 2>&1
  exec > >(tee test.log)
  echo "$@"
  echo something
}

alirootInfo()
{
  (
  umask 0002
  # save aliroot svn info
  [[ -z "$ALICE_ROOT" ]] && exit 1
  prevdir=$PWD
  cd $ALICE_ROOT
  echo "\$ALICE_ROOT=$ALICE_ROOT"
  echo
  svn info 2>/dev/null
  echo ""
  echo ""
  svn diff 2>/dev/null
  cd $prevdir
  )
}

MakeListOfQAresults()
{
  (
  umask 0002
  if [[ $# -eq 0 ]]; then
    echo "make an input file list for the qa script"
    echo "Usage: $0 path filename"
    return 0
  fi

  path=$1
  fileName=$2
  
  while read entry; do
    runNumber=$(guessRunNumber $entry)
    echo "$runNumber 0 1 10 $entry"
  done < <(find $path -name $fileName)
  )
}

setYear()
{
  #set the year
  #  $1 - year to be set
  #  $2 - where to set the year
  year1=$(guessYear $1)
  year2=$(guessYear $2)
  local path=$2
  [[ $year1 -ne $year2 && -n $year2 ]] && path=${2/\/$year2\//\/$year1\/}
  echo $path
}

guessPeriod()
{
  #guess the period from the path, pick the rightmost one
  local IFS="/"
  declare -a path=( $1 )
  local dirDepth=${#path[*]}
  for ((x=${dirDepth}-1;x>=0;x--)); do
    local field=${path[${x}]}
    [[ ${field} =~ ^LHC[0-9][0-9][a-z]$ ]] && period=${field#000} && break
  done
  echo $period
}

guessYear()
{
  #guess the year from the path, pick the rightmost one
  local IFS="/"
  declare -a path=( $1 )
  local dirDepth=${#path[*]}
  for ((x=${dirDepth}-1;x>=0;x--)); do
    local field=${path[${x}]}
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field#000} && break
  done
  echo $year
}

guessRunNumber()
{
  #guess the run number from the path, pick the rightmost one
  local IFS="/"
  declare -a path=( $1 )
  local dirDepth=${#path[*]}
  for ((x=${dirDepth}-1;x>=0;x--)); do
    local field=${path[${x}]}
    [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && runNumber=${field#000} && break
  done
  echo $runNumber
}

summarizeLogs()
{
  #print a summary of logs
  logFiles=(
            "*.log"
  )

  errorConditions=(
                    "There was a crash"
                    "floating"
                    "error while loading shared libraries"
                    "std::bad_alloc"
                    "s_err_syswatch_"
                    "Thread [0-9]* (Thread"
  )
  logstatus=0
  for log in ${logFiles[*]}; do
    finallog=${outputDir%/}/${log}
    [[ ! -f $log ]] && continue
    errorSummary=""
    for ((i=0; i<${#errorConditions[@]};i++)); do
      local tmp=$(grep -m1 -e "${errorConditions[$i]}" $log)
      [[ -n $tmp ]] && tmp+=" : "
      errorSummary+=$tmp
    done
    if [[ -z $errorSummary ]]; then 
      #in pretend mode randomly report an error in rec.log some cases
      if [[ -n $pretend && "$log" == "rec.log" ]]; then
        [[ $(( $RANDOM%2 )) -ge 1 ]] && echo "$finallog BAD random error" || echo "$finallog OK"
      else
        echo "$finallog OK"
      fi
    else 
      local logstatus=1
      echo "$finallog BAD $errorSummary"
    fi
  done
  return $logstatus
}

spitOutLocalOCDBaccessConfig()
{
  umask 0002
  find $1 -name "*root" | \
  while read line
  do 
    local tmp=${line#$1}
    echo ${tmp%/*} | \
    awk -v ocdb=$1 '{print "  man->SetSpecificStorage(\""$1"\",\"local://"ocdb"\");"}'
  done
}

goConfigureOCDBforCPass1()
{
  umask 0002
  # make a script that sets the specific storages form all the root files produced by CPass0
  # in this case the second argument is taken to be the path to the produced OCDB for specific storage
  local localOCDBpathCPass0=$1
  local outputDir=$2
  [[ -d $outputDir ]] && cd $outputDir

  if [[ -f $localOCDBpathCPass0 && $localOCDBpathCPass0 =~ \.tgz$ ]]; then
    tar xzf $localOCDBpathCPass0
    local localOCDBpathCPass0="./OCDB"
  fi

  echo
  echo creating the specific storage script
  echo   localOCDBaccessConfig.C
  echo   based on OCDB: $fileName
  echo

  local tempLocalOCDB=""
  if [[ -f localOCDBaccessConfig.C ]]; then
    tempLocalOCDB=$(mktemp)
    echo "egrep "SetSpecificStorage" localOCDBaccessConfig.C > $tempLocalOCDB"
    egrep "SetSpecificStorage" localOCDBaccessConfig.C > $tempLocalOCDB
  fi

  echo "localOCDBaccessConfig()"                               >  localOCDBaccessConfig.C
  echo "{"                                                     >> localOCDBaccessConfig.C
  echo "  AliCDBManager* man = AliCDBManager::Instance();"     >> localOCDBaccessConfig.C
  spitOutLocalOCDBaccessConfig $localOCDBpathCPass0|sort|uniq  >> localOCDBaccessConfig.C
  [[ -f "$tempLocalOCDB" ]] && cat $tempLocalOCDB              >> localOCDBaccessConfig.C
  echo "}"                                                     >> localOCDBaccessConfig.C

  [[ -f "$tempLocalOCDB" ]] && rm -f $tempLocalOCDB

  if ! grep SetSpecificStorage localOCDBaccessConfig.C; then 
    echo
    echo "!!!!!!! CPass0 produced no OCDB entries"
    return 1
  fi
}

goSubmit()
{
  inputList=$1
  productionID=$2
  configFile="benchmark.config"
  [[ -n "$3" ]] && configFile=$3
  configFile=$(readlink -f $configFile)
  [[ -n $4 ]] && runNumber=$4

  #redirect all output to submit.log
  echo "redirecting all output to ${PWD}/submit_${productionID//"/"/_}.log"
  exec 7>&1
  exec 1>submit_${productionID//"/"/_}.log 2>&1

  umask 0002
  echo $0" submit $*"
  if [[ -z "$inputList" || -z "$productionID" ]]
  then
    echo
    echo " Usage: $0 submit inputList productionID [configFile=benchmark.config]"
    echo
    exit
  fi

  date=`date +%Y-%m-%d_%H%M%S`

  # check if config file is there
  if [ ! -f $configFile ]; then
    echo "ERROR! Config File '$configFile' not found" >&2
    exit
  else
    echo "Using Config File: '$configFile'"
  fi

  # source the config file
  # and print the configuration
  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} && -z ${alirootEnv} ]] && source $setupAliROOTenvInCurrentShell

  self=$(readlink -f "$0")
  configPath=`dirname $self`
  #exporting makes it available in the environment on the nodes - makes the summary output go there
  export commonOutputPath=${baseOutputDirectory}/${productionID}

  #[[ -z ${ALIROOT_RELEASE} ]] && ALIROOT_RELEASE=${ALICE_ROOT//"/"/"_"}
  #ALIROOT_RELEASE_mod=${ALIROOT_RELEASE//"/"/"_"}
  #ALIROOT_BASEDIR_mod=${ALIROOT_BASEDIR//"/"/"_"}
  #productionID=${ALIROOT_BASEDIR_mod}_${ALIROOT_RELEASE_mod}/${productionID}

  #convert to absolut pathnames
  inputList=$(readlink -f "$inputList")
  #make list of runs
  if [[ -z $runNumber ]]; then
    listOfRuns=($(while read x; do guessRunNumber $x; done < ${inputList} | sort | uniq))
  else
    listOfRuns=$runNumber
  fi

  setupAliROOTenvInCurrentShell=$(readlink -f "$setupAliROOTenvInCurrentShell")

  echo ""
  echo "### BEGIN CONFIGURATION ###"
  echo ""
  echo "GENERAL:"
  echo ""
  echo "    productionID:    $productionID"
  echo "    batchCommand:    $batchCommand"
  echo "    setupAliROOTenvInCurrentShell:      $setupAliROOTenvInCurrentShell"
  echo "    ALICE_ROOT:      $ALICE_ROOT"
  echo "    ALIROOT_RELEASE: $ALICE_RELEASE"
  echo "    inputList:       $inputList"
  echo "    configPath:      $configPath"
  echo "    commonOutputPath:      $commonOutputPath"
  echo "    defaultOCDB:     $defaultOCDB"
  echo "      autoOCDB: $autoOCDB"
  echo "    recoTriggerOptions:   $recoTriggerOptions"
  echo "    runs:"
  echo "      ${listOfRuns[*]}"
  echo ""
  echo "THE TRAIN WILL RUN:"

  if [ $runCPass0reco -eq 1 ]; then
    echo "    Pass0 - Recontruction"
  fi

  if [ $runCPass0MergeMakeOCDB -eq 1 ]; then
    echo "    Pass0 - merging and OCDB export"
  fi

  if [ $runCPass1reco -eq 1 ]; then
    echo "    Pass1 - Recontruction"
  fi
  if [ $runCPass1MergeMakeOCDB -eq 1 ]; then
    echo "    Pass1 - merging and OCDB export"
  fi

  echo ""
  echo "LIMITS:"
  echo "    max. Events/Chunk:   $nEvents"
  echo "    max. Number of Chunks per Run:     $nMaxChunks"
  echo ""
  echo "### END CONFIGURATION ###"
  echo ""


  # check if input file is there
  if [ ! -f $inputList ]; then
    echo "ERROR! Input List '$inputList' not found" >&2
    exit
  fi

  # define jobid (for dependent jobs)
  JOBID1=p0_${productionID//"/"/_}_${date}
  JOBID1wait=w0_${productionID//"/"/_}_${date}
  JOBID2=m0_${productionID//"/"/_}_${date}
  JOBID2wait=wm0_${productionID//"/"/_}_${date}
  JOBID3=op0_${productionID//"/"/_}_${date}
  JOBID3wait=wop0_${productionID//"/"/_}_${date}
  JOBID4=p1_${productionID//"/"/_}_${date}
  JOBID4wait=w1_${productionID//"/"/_}_${date}
  JOBID5=m1_${productionID//"/"/_}_${date}
  JOBID5wait=wm1_${productionID//"/"/_}_${date}
  JOBID6=s1_${productionID//"/"/_}_${date}
  JOBID6wait=ws1_${productionID//"/"/_}_${date}
  JOBID7=QA_${productionID//"/"/_}_${date}
  LASTJOB=""

  #for each run we submit one jobarray:
  numberOfSubmittedCPass1MergingJobs=0
  for runNumber in ${listOfRuns[*]}; do
    [[ -z $runNumber ]] && continue
    [[ ! ${runNumber} =~ ^[0-9]*[0-9]$ ]] && continue
    oneInputFile=$(egrep -m1 "$runNumber\/" ${inputList})
    currentDefaultOCDB=$defaultOCDB
    [[ $autoOCDB -ne 0 ]] && currentDefaultOCDB=$(setYear $oneInputFile $defaultOCDB)
    echo "submitting run $runNumber with OCDB $currentDefaultOCDB"

    ################################################################################
    ################################################################################
    # run the CPass0 if requested

    if [ $runCPass0reco -eq 1 ]; then

      echo
      echo "starting CPass0... for run $runNumber"
      echo

      # create directory and copy all files that are needed
      targetDirectory="${commonOutputPath}/000${runNumber}/CPass0"
      mkdir -p $targetDirectory
      mkdir -p $targetDirectory/logs
      
      localInputList=$targetDirectory/${inputList##*/}
      theScript="$targetDirectory/${self##*/}"
      cp -f $self $theScript
      chmod u+x $theScript

      rm -f $localInputList
      egrep "\/000$runNumber\/" $inputList >> $localInputList
      cp -f $configFile $targetDirectory
      [[ -f $configPath/runCPass0.sh ]] && cp -f $configPath/runCPass0.sh $targetDirectory && echo "## using local runCPass0.sh"
      [[ -f $configPath/recCPass0.C ]] && cp -f $configPath/recCPass0.C $targetDirectory && echo "## using local recCPass0.C"
      [[ -f $configPath/runCalibTrain.C ]] && cp -f $configPath/runCalibTrain.C $targetDirectory && echo "## using local runCalibTrain.C"
      [[ -f $configPath/localOCDBaccessConfig.C ]] && cp -f $configPath/localOCDBaccessConfig.C $targetDirectory && echo "## using local localOCDBaccessConfig.C"

      [[ -f $configPath/CPass0/runCPass0.sh ]] && cp -f $configPath/CPass0/runCPass0.sh $targetDirectory && echo "## using local runCPass0.sh"
      [[ -f $configPath/CPass0/recCPass0.C ]] && cp -f $configPath/CPass0/recCPass0.C $targetDirectory && echo "## using local recCPass0.C"
      [[ -f $configPath/CPass0/runCalibTrain.C ]] && cp -f $configPath/CPass0/runCalibTrain.C $targetDirectory && echo "## using local runCalibTrain.C"
      [[ -f $configPath/CPass0/localOCDBaccessConfig.C ]] && cp -f $configPath/CPass0/localOCDBaccessConfig.C $targetDirectory && echo "## using local localOCDBaccessConfig.C"

      echo "... files copied."
      echo

      # limit nFiles to nMaxChunks
      nFiles=`wc -l < $localInputList`
      [[ $nFiles -eq 0 ]] && echo "list contains ZERO files! exiting..." && exit 1
      echo "raw files in list:    $nFiles"
      if [[ $nMaxChunks -gt 0 && $nMaxChunks -le $nFiles ]]; then
        nFiles=$nMaxChunks
      fi
      echo "raw files to process: $nFiles"
      [[ -z "$percentProcessedFilesToContinue" ]] && percentProcessedFilesToContinue=100
      if [[ $percentProcessedFilesToContinue -eq 100 ]]; then
        nFilesToWaitFor=$nFiles
      else
        nFilesToWaitFor=$(( $nFiles-$nFiles/(100/(100-$percentProcessedFilesToContinue)) ))
      fi
      echo "requested success rate is $percentProcessedFilesToContinue%"
      echo "merging will start after $nFilesToWaitFor jobs are done"

      submit $JOBID1 1 $nFiles "" "$theScript" "CPass0 $targetDirectory $localInputList $nEvents $currentDefaultOCDB $configFile $runNumber"

      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "$JOBID1wait" 1 1 "" "$theScript" "WaitForOutput ${commonOutputPath} 'cpass0.job*.run$runNumber.done' $nFilesToWaitFor $maxSecondsToWait '-maxdepth 1'"
      LASTJOB=$JOBID1wait

    fi #end running CPass0
    ################################################################################


    ################################################################################
    # submit merging of CPass0, depends on the reconstruction

    if [ $runCPass0MergeMakeOCDB -eq 1 ]; then

      echo
      echo "submit CPass0 merging for run $runNumber"
      echo

      targetDirectory="${commonOutputPath}/000${runNumber}/CPass0"
      mkdir -p $targetDirectory

      # copy the scripts
      cp -f $self $targetDirectory
      [[ -f $configPath/mergeMakeOCDB.sh ]] && cp -f $configPath/mergeMakeOCDB.sh $targetDirectory && echo "## using local mergeMakeOCDB.sh"
      [[ -f $configPath/merge.C ]] && cp -f $configPath/merge.C $targetDirectory && echo "## using local merge.C"
      [[ -f $configPath/mergeMakeOCDB.byComponent.sh ]] && cp -f $configPath/mergeMakeOCDB.byComponent.sh $targetDirectory && echo "## using local mergeMakeOCDB.byComponent.sh"
      [[ -f $configPath/mergeByComponent.C ]] && cp -f $configPath/mergeByComponent.C $targetDirectory && echo "## using local mergeByComponent.C"
      [[ -f $configPath/makeOCDB.C ]] && cp -f $configPath/makeOCDB.C $targetDirectory && echo "## using local makeOCDB.C"

      theScript="$targetDirectory/${self##*/}"
      chmod u+x $theScript

      echo submit $JOBID2 1 1 "$LASTJOB" "$theScript" "MergeCPass0 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      submit $JOBID2 1 1 "$LASTJOB" "$theScript" "MergeCPass0 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      LASTJOB=$JOBID2

      cd $configPath
      echo
    fi
    # end of merging CPass0
    ################################################################################

    ################################################################################
    ################################################################################
    # run the CPass1 if requested

    if [ $runCPass1reco -eq 1 ]; then

      targetDirectory="${commonOutputPath}/000${runNumber}/CPass1"

      # safety feature: if we are re-running for any reason we want to delete the previous output first.
      [[ -d $targetDirectory ]] && rm -rf $targetDirectory/* && echo "removed old output at $targetDirectory/*"

      ################################################################################
      # for the CPass1, based on the OCDB entries produced, we need to create the script 
      # to set the specific storages to the output of CPass0
      # submit a job that will execute after merging/OCDB export that will do it.

      mkdir -p $targetDirectory
      mkdir -p $targetDirectory/logs
      cp -f $self $targetDirectory
      theScript="$targetDirectory/${self##*/}"

      echo
      echo submitting the OCDB specific storage config making script for run $runNumber
      echo

      [[ -f $configPath/localOCDBaccessConfig.C ]] && cp -f $configPath/localOCDBaccessConfig.C $targetDirectory && echo "## using local localOCDBaccessConfig.C"

      submit $JOBID3 1 1 "$LASTJOB" "$theScript" "ConfOCDB ${commonOutputPath}/000${runNumber}/CPass0/OCDB ${targetDirectory} "
      LASTJOB=$JOBID3
      ################################################################################

      echo
      echo "starting CPass1... for run $runNumber"
      echo

      # create directory and copy all files that are needed
      mkdir -p $targetDirectory
      mkdir -p $targetDirectory/logs
      
      localInputList=$targetDirectory/${inputList##*/}
      theScript="$targetDirectory/${self##*/}"
      cp -f $self $theScript
      chmod u+x $theScript

      rm -f $localInputList
      egrep "\/000$runNumber\/" $inputList >> $localInputList
      cp -f $configFile $targetDirectory
      [[ -f $configPath/runCPass1.sh ]] && cp -f $configPath/runCPass1.sh $targetDirectory && echo "## using local runCPass1.sh"
      [[ -f $configPath/recCPass1.C ]] && cp -f $configPath/recCPass1.C $targetDirectory && echo "## using local recCPass1.C"
      [[ -f $configPath/recCPass1_OuterDet.C ]] && cp -f $configPath/recCPass1_OuterDet.C $targetDirectory && echo "## using local recCPass1_OuterDet.C"
      [[ -f $configPath/runCalibTrain.C ]] && cp -f $configPath/runCalibTrain.C $targetDirectory && echo "## using local runCalibTrain.C"
      [[ -f $configPath/QAtrain.C ]] && cp -f $configPath/QAtrain.C $targetDirectory && echo "## using local QAtrain.C"
      [[ -f $configPath/QAtrain_duo.C ]] && cp -f $configPath/QAtrain_duo.C $targetDirectory && echo "## using local QAtrain_duo.C"

      [[ -f $configPath/CPass1/runCPass1.sh ]] && cp -f $configPath/CPass1/runCPass1.sh $targetDirectory && echo "## using local runCPass1.sh"
      [[ -f $configPath/CPass1/recCPass1.C ]] && cp -f $configPath/CPass1/recCPass1.C $targetDirectory && echo "## using local recCPass1.C"
      [[ -f $configPath/CPass1/recCPass1_OuterDet.C ]] && cp -f $configPath/CPass1/recCPass1_OuterDet.C $targetDirectory && echo "## using local recCPass1_OuterDet.C"
      [[ -f $configPath/CPass1/runCalibTrain.C ]] && cp -f $configPath/CPass1/runCalibTrain.C $targetDirectory && echo "## using local runCalibTrain.C"
      [[ -f $configPath/CPass1/QAtrain.C ]] && cp -f $configPath/CPass1/QAtrain.C $targetDirectory && echo "## using local QAtrain.C"
      [[ -f $configPath/CPass1/QAtrain_duo.C ]] && cp -f $configPath/CPass1/QAtrain_duo.C $targetDirectory && echo "## using local QAtrain_duo.C"

      echo "... files copied."
      echo

      # limit nFiles to nMaxChunks
      nFiles=`wc -l < $localInputList`
      [[ $nFiles -eq 0 ]] && echo "list contains ZERO files! exiting..." && exit 1
      echo "raw files in list:    $nFiles"
      if [[ $nMaxChunks -gt 0 && $nMaxChunks -le $nFiles ]]; then
        nFiles=$nMaxChunks
      fi
      echo "raw files to process: $nFiles"
      [[ -z "$percentProcessedFilesToContinue" ]] && percentProcessedFilesToContinue=100
      if [[ $percentProcessedFilesToContinue -eq 100 ]]; then
        nFilesToWaitFor=$nFiles
      else
        nFilesToWaitFor=$(( $nFiles-$nFiles/(100/(100-$percentProcessedFilesToContinue)) ))
      fi
      echo "requested success rate is $percentProcessedFilesToContinue%"
      echo "merging will start after $nFilesToWaitFor jobs are done"

      submit $JOBID4 1 $nFiles "$LASTJOB" "$theScript" "CPass1 $targetDirectory $localInputList $nEvents $currentDefaultOCDB $configFile $runNumber"

      ################################################################################
      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "$JOBID4wait" 1 1 "$LASTJOB" "$theScript" "WaitForOutput ${commonOutputPath} 'cpass1.job*.run${runNumber}.done' $nFilesToWaitFor $maxSecondsToWait '-maxdepth 1'"
      LASTJOB=$JOBID4wait
      ################################################################################

      echo
    fi #end running CPass1

    ################################################################################
    # submit merging of CPass1, depends on the reconstruction
    if [ $runCPass1MergeMakeOCDB -eq 1 ]; then

      echo
      echo "submit CPass1 merging for run $runNumber"
      echo

      targetDirectory="${commonOutputPath}/000${runNumber}/CPass1"
      mkdir -p $targetDirectory

      # copy the scripts
      cp -f $self $targetDirectory
      [[ -f $configPath/mergeMakeOCDB.sh ]] && cp -f $configPath/mergeMakeOCDB.sh $targetDirectory && echo "## using local mergeMakeOCDB.sh"
      [[ -f $configPath/merge.C ]] && cp -f $configPath/merge.C $targetDirectory && echo "## using local merge.C"
      [[ -f $configPath/mergeMakeOCDB.byComponent.sh ]] && cp -f $configPath/mergeMakeOCDB.byComponent.sh $targetDirectory && echo "## using local mergeMakeOCDB.byComponent.sh"
      [[ -f $configPath/mergeByComponent.C ]] && cp -f $configPath/mergeByComponent.C $targetDirectory && echo "## using local mergeByComponent.C"
      [[ -f $configPath/makeOCDB.C ]] && cp -f $configPath/makeOCDB.C $targetDirectory && echo "## using local makeOCDB.C"

      theScript="$targetDirectory/${self##*/}"
      chmod u+x $theScript

      echo submit "$JOBID5" 1 1 "$LASTJOB" "$theScript" "MergeCPass1 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      submit "$JOBID5" 1 1 "$LASTJOB" "$theScript" "MergeCPass1 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      ((numberOfSubmittedCPass1MergingJobs++))
      LASTJOB=$JOBID5
      echo
    fi

  done

  ################################################################################
  ################################################################################
  [[ -z $runMakeSummary ]] && runMakeSummary=0
  if [ $runMakeSummary -eq 1 ]; then
    echo
    echo "submit make a summary"
    echo

    targetDirectory="${commonOutputPath}"
    mkdir -p $targetDirectory
    mkdir -p $targetDirectory/logs
    [[ ! -f $targetDirectory/${self##*/} ]] && cp -f $self $targetDirectory
    theScript="$targetDirectory/${self##*/}"
    submit "$JOBID6" 1 1 "$LASTJOB" "$theScript" "MakeSummary $targetDirectory $configFile"
  fi
  ################################################################################
  
  ################################################################################
  ################################################################################
  if [ $runMakeQAplots -eq 1 ]; then
    echo
    echo "submit make QA plots"
    echo

    ## submit a monitoring job that will wait for a specified number of files to appear somewhere
    targetDirectory="${commonOutputPath}"
    [[ ! -f $targetDirectory/${self##*/} ]] && cp -f $self $targetDirectory
    theScript="$targetDirectory/${self##*/}"
    echo submit "$JOBID5wait 1 1 $theScript WaitForOutput ${commonOutputPath} 'merge.cpass1.run*.done' $numberOfSubmittedCPass1MergingJobs $maxSecondsToWait '-maxdepth 1'"
    submit "$JOBID5wait" 1 1 "" "$theScript" "WaitForOutput ${commonOutputPath} 'merge.cpass1.run*.done' $numberOfSubmittedCPass1MergingJobs $maxSecondsToWait '-maxdepth 1'"
    LASTJOB=$JOBID5wait
  
    targetDirectory="${commonOutputPath}/QAplots/"
    [[ -d $targetDirectory ]] && rm -rf $targetDirectory/*
    mkdir -p $targetDirectory
    mkdir -p $targetDirectory/logs
    qaFilesDirectory="${commonOutputPath}"
    [[ ! -f $targetDirectory/${self##*/} ]] && cp -f $self $targetDirectory
    theScript="$targetDirectory/${self##*/}"

    echo submit "$JOBID7" 1 1 "$LASTJOB" "$theScript" "CreateQAplots QAplots CPass1 $targetDirectory $qaFilesDirectory $qaPlotsScript $configFile"
    submit "$JOBID7" 1 1 "$LASTJOB" "$theScript" "CreateQAplots QAplots CPass1 $targetDirectory $qaFilesDirectory $qaPlotsScript $configFile"
    LASTJOB=$JOBID7
  fi

  #restore stdout
  exec 1>&7 7>&-
  echo "jobs submitted."
}

main "$@"
