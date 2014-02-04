#!/bin/bash
# this script runs the CPass0/CPass1 train
# produced OCDB updates are local

main()
{
  #run in proper mode depending on the selection
  runMode=$1
  #if [[ $# -gt 0 ]]; then
  #  echo "# $0 $*"
  #fi
  umask 0002
  shift
  case $runMode in
    "CPass0") goCPass0 "$@";;
    "CPass1") goCPass1 "$@";;
    "ConfOCDB") goConfigureOCDBforCPass1 "$@";;
    "MergeCPass0") goMergeCPass0 "$@";;
    "MergeCPass1") goMergeCPass1 "$@";;
    "MakeFilteredTrees") goMakeFilteredTrees "$@";;
    "makeSummary") goMakeSummary "$@";;
    "run") goSubmitMakeflow "$@";;
    "submit") goSubmitBatch "$@";;
    "test") goTest "$@";;
    "generateMakeflow") goGenerateMakeflow "$@";;
    "printValues") goPrintValues "$@";;
    "CreateQAplots") goCreateQAplots "$@";;
    "WaitForOutput") goWaitForOutput "$@";;
    "merge") goMerge "$@";;
    #and the default: unless sourced print some info
    *) if [[ ! "$0" =~ "bash" ]]; then
         echo "uses makeflow:"
         echo " $0 \"run\" inputList productionID [configFile=benchmark.config] [runNumber]"
         echo "uses a batch system (SGE):"
         echo " $0 \"submit\" inputList productionID [configFile=benchmark.config] [runNumber]"
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

  #the contents of this is stored in the tree and used later (e.g. AliAnalysisTaskPIDResponse)!
  #at the QA stage the pass number is guessed from the path stored here.
  #The Format is:
  #Packages= ;OutputDir= ;LPMPass= ;TriggerAlias= ;LPMRunNumber= ;LPMProductionType= ;LPMInteractionType= ;LPMProductionTag= ;LPMAnchorRun= ;LPMAnchorProduction= ;LPMAnchorYear= 
  export PRODUCTION_METADATA="OutputDir=cpass0"

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
  outputDir=${targetDirectory}/${jobindex}_${chunkName}
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
    touch AliESDs.root
    touch AliESDfriends.root
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
  cp -p --recursive $runpath/* $outputDir
  echo
  
  #validate CPass0
  cd ${outputDir}
  touch ${calibDoneFile}
  [[ -f AliESDfriends_v1.root ]] && echo "calibfile ${outputDir}/AliESDfriends_v1.root" > ${calibDoneFile}
  [[ -f AliESDs.root ]] && echo "esd ${outputDir}/AliESDs.root" >> ${calibDoneFile}
  echo "dir ${outputDir}" >> ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}

  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
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

  #the contents of this is stored in the tree and used later (e.g. AliAnalysisTaskPIDResponse)!
  #at the QA stage the pass number is guessed from the path stored here.
  #The Format is:
  #Packages= ;OutputDir= ;LPMPass= ;TriggerAlias= ;LPMRunNumber= ;LPMProductionType= ;LPMInteractionType= ;LPMProductionTag= ;LPMAnchorRun= ;LPMAnchorProduction= ;LPMAnchorYear= 
  export PRODUCTION_METADATA="OutputDir=cpass1"

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
  outputDir=${targetDirectory}/${jobindex}_${chunkName}
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
               "$trustedQAtrainMacro"
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

  if [[ ! $(find ./OCDB -name "*.root") ]]; then
    echo "cpass0 produced no calibration! exiting..."
    exit 1
  fi

  #run CPass1
  chmod u+x runCPass1.sh
  echo "$runpath/runCPass1.sh $infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
  if [[ -n $pretend ]]; then
    touch AliESDs_Barrel.root
    touch AliESDfriends_Barrel.root
    touch AliESDfriends_v1.root
    touch QAresults_barrel.root
    touch EventStat_temp_barrel.root
    touch AODtpITS.root
    touch AliESDs_Outer.root
    touch AliESDfriends_Outer.root
    touch QAresults_outer.root
    touch EventStat_temp_outer.root
    touch rec.log
    touch calib.log
    touch qa.log
  else
    ./runCPass1.sh "$infile" "$nEvents" "$runNumber" "$ocdbPath" "$recoTriggerOptions"
  fi

  ##handle possible crashes in QA (happens often in trunk)
  ##rerun QA with a trusted aliroot version
  #if [[ $(validateLog qa_barrel.log) ]]; then
  #  echo "qa_barrel.log not validated!"
  #fi
  #if [[ ! -f QAresults_barrel.root && -f ${setupTrustedAliROOTenvInCurrentShell} || $(validateLog qa_barrel.log) ]]; then
  #  echo "WARNING: using trusted QA aliroot $ALICE_ROOT"
  #  source ${setupTrustedAliROOTenvInCurrentShell}
  #  cd Barrel
  #  rm QAresults_barrel.root
  #  rm EventStat_temp_barrel.root
  #  rm AODtpITS.root
  #  [[ ! -f AliESDs.root ]] && ln -s ../AliESDs_Barrel.root AliESDs.root
  #  [[ ! -f AliESDfriends.root ]] && ln -s ../AliESDfriends_Barrel.root AliESDfriends.root
  #  if [[ -n $trustedQAtrainMacro ]]; then
  #    eval "cp $trustedQAtrainMacro QAtrain_duo_trusted.C"
  #  fi
  #  echo executing aliroot -b -q "QAtrain_duo_trusted.C(\"_barrel\",$runNumber,\"wn.xml\",0,\"$ocdbPath\")"
  #  time aliroot -b -q "QAtrain_duo.C(\"_barrel\",$runNumber,\"wn.xml\",0,\"$ocdbPath\")" &> ../qa_barrel_trusted.log
  #  cd ../
  #fi
  [[ ! -f AliESDs_Barrel.root && -f Barrel/AliESDs.root ]] && mv Barrel/AliESDs.root AliESDs_Barrel.root
  [[ ! -f AliESDfriends_Barrel.root && -f Barrel/AliESDfriends.root ]] && mv Barrel/AliESDfriends.root AliESDfriends_Barrel.root
  [[ ! -f AliESDfriends_v1.root && -f Barrel/AliESDfriends_v1.root ]] && mv Barrel/AliESDfriends_v1.root .
  [[ ! -f QAresults_barrel.root && -f Barrel/QAresults_barrel.root ]] && mv Barrel/QAresults_barrel.root .
  [[ ! -f AliESDs_Outer.root && -f OuterDet/AliESDs.root ]] && mv OuterDet/AliESDs.root AliESDs_Outer.root
  [[ ! -f AliESDfriends_Outer.root && -f OuterDet/AliESDfriends.root ]] && mv OuterDet/AliESDfriends.root AliESDfriends_Outer.root
  [[ ! -f QAresults_outer.root && -f OuterDet/QAresults_outer.root ]] && mv OuterDet/QAresults_outer.root .

  #move stuff to final destination
  echo "this directory ($PWD) contents:"
  ls
  echo rm -f ./$chunkName
  rm -f ./$chunkName
  echo "cp --recursive ${runpath}/* ${outputDir}"
  cp -p --recursive ${runpath}/* ${outputDir}
  echo

  #validate CPass1
  cd ${outputDir}
  touch ${calibDoneFile}
  [[ -f AliESDs_Barrel.root ]] && echo "esd ${outputDir}/AliESDs_Barrel.root" > ${calibDoneFile}
  [[ -f AliESDfriends_v1.root ]] && echo "calibfile ${outputDir}/AliESDfriends_v1.root" >> ${calibDoneFile}
  [[ -f QAresults_Barrel.root ]] && echo "qafile ${outputDir}/QAresults_Barrel.root" >> ${calibDoneFile}
  [[ -f QAresults_Outer.root ]] && echo "qafile ${outputDir}/QAresults_Outer.root" >> ${calibDoneFile}
  echo "dir ${outputDir}" >> ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}
  
  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
}


goMergeCPass0()
{
  umask 0002
  #
  # find the output files and merge them
  #

  outputDir=$1
  ocdbStorage=$2
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
  [[ -z $dontRedirectStdOutToLog ]] && exec &> $logOutputDir/mergeMakeOCDB.log
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
  echo ocdbStorage=$ocdbStorage
  echo calibrationFilesToMerge=$calibrationFilesToMerge
  echo calibrationOutputFileName=$calibrationOutputFileName
  echo mergingScript=$mergingScript
  echo commonOutputPath=$commonOutputPath
  echo runpath=$runpath
  
  # copy files in case they are not already there
  filesMergeCPass0=(
                    "$commonOutputPath/${calibrationFilesToMerge}"
                    "$commonOutputPath/OCDB.root"
                    "$commonOutputPath/localOCDBaccessConfig.C"
                    "$commonOutputPath/cpass0.localOCDB.${runNumber}.tgz"
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
  echo "PWD"
  ls -lh
  echo "PWD/.."
  ls -lh ../


  #merge calibration
  chmod u+x $mergingScript  
  mkdir -p ./OCDB
  if [[ -z ${calibrationFilesToMergeExternal} ]]; then
    echo "find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge"
    find $outputDir -name $calibrationOutputFileName > $calibrationFilesToMerge
  fi
  
  echo "$mergingScript $calibrationFilesToMerge ${runNumber} local://./OCDB $ocdbStorage"
  if [[ -n $pretend ]]; then
    touch CalibObjects.root
    touch ocdb.log
    touch merge.log
    mkdir -p ./OCDB/someDetector/
    mkdir -p ./OCDB/otherDetector/
    touch ./OCDB/someDetector/someCalibObject_0-999999_cpass0.root
    touch ./OCDB/otherDetector/otherCalibObject_0-999999_cpass0.root
  else
    ./$mergingScript $calibrationFilesToMerge ${runNumber} "local://./OCDB" $ocdbStorage
  fi

  ### produce the output
  #tar the produced OCDB for reuse
  tar czf $commonOutputPath/cpass0.localOCDB.${runNumber}.tgz ./OCDB

  ls -ltrh

  #copy all to output dir
  cp -p --recursive ${runpath}/* $outputDir
  
  #validate merging cpass0
  cd ${outputDir}
  calibDoneFile="${commonOutputPath}/merge.cpass0.run${runNumber}.done"
  touch ${calibDoneFile}
  [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" > ${calibDoneFile}
  summarizeLogs >> ${calibDoneFile}

  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
}

goMergeCPass1()
{
  umask 0002
  #
  # find the output files and merge them
  #

  outputDir=$1
  ocdbStorage=$2
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
  [[ -z $dontRedirectStdOutToLog ]] && exec &> $logOutputDir/mergeMakeOCDB.log
  echo "$0 $*"

  calibrationFilesToMerge=$calibrationFilesToMergeExternal
  [[ -z $calibrationFilesToMerge ]] && calibrationFilesToMerge="calibrationFilesToMerge.list"
  calibrationOutputFileName="AliESDfriends_v1.root"
  mergingScript="mergeMakeOCDB.byComponent.sh"
  qaFilesToMerge=$qaFilesToMergeExternal
  #important to have the string "Stage.txt" in the filename to trigger the merging
  #it has to be a list of directories containing the files
  [[ -z $qaFilesToMerge ]] && qaFilesToMerge="qaFilesForMergingStage.txt.list"
  qaOutputFileName="QAresults*.root"
  qaMergedOutputFileName="QAresults_merged.root"

  echo goMergeCPass1 SETUP:
  echo runNumber=$runNumber
  echo outputDir=$outputDir
  echo ocdbStorage=$ocdbStorage
  echo calibrationFilesToMerge=filesToMerge.list
  echo calibrationOutputFileName=$calibrationOutputFileName
  echo mergingScript=$mergingScript
  
  # copy files in case they are not already there
  filesMergeCPass1=(
                    "$commonOutputPath/${calibrationFilesToMerge}"
                    "$commonOutputPath/${qaFilesToMerge}"
                    "$commonOutputPath/OCDB.root"
                    "$commonOutputPath/localOCDBaccessConfig.C"
                    "$commonOutputPath/cpass0.localOCDB.${runNumber}.tgz"
                    "$commonOutputPath/QAtrain_duo.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.byComponent.sh"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeByComponent.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/makeOCDB.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/merge.C"
                    "$ALICE_ROOT/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.sh"
                    "$ALICE_ROOT/ANALYSIS/macros/QAtrain_duo.C"
  )
  for file in ${filesMergeCPass1[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done

  #configure local OCDB storage from CPass0 (creates the localOCDBaccessConfig.C script)
  if [[ -f cpass0.localOCDB.${runNumber}.tgz ]]; then
    echo goConfigureOCDBforCPass1 "cpass0.localOCDB.${runNumber}.tgz"
    goConfigureOCDBforCPass1 "cpass0.localOCDB.${runNumber}.tgz"
  else
    echo "WARNING: file cpass0.localOCDB.${runNumber}.tgz not found!"
  fi

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
  
  echo "$mergingScript $calibrationFilesToMerge ${runNumber} local://./OCDB $ocdbStorage"
  if [[ -n $pretend ]]; then
    touch CalibObjects.root
    touch ocdb.log
    touch merge.log
  else
    ./$mergingScript $calibrationFilesToMerge ${runNumber} "local://./OCDB" $ocdbStorage
  fi

  tar czf localCPass1_${runNumber}.tgz ./OCDB

  #merge QA
  [[ -n ${AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF
  [[ -n ${AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF

  if [[ -z $qaFilesToMergeExternal ]]; then
    #find the files, but only store the directories (QAtrain_duo.C requires this)
    echo "find $outputDir -name $qaOutputFileName | while read x; do echo ${x%/*}; done | sort | uniq > $qaFilesToMerge"
    find $outputDir -name $qaOutputFileName | while read x; do echo ${x%/*}; done | sort | uniq > $qaFilesToMerge
  fi
  
  #do the merge
  #echo aliroot -l -b -q "merge.C(\"$qaFilesToMerge\",\"\",kFALSE,\"$qaMergedOutputFileName\")"
  echo aliroot -b -q "QAtrain_duo.C(\"_barrel\",${runNumber},\"$qaFilesToMerge\",1,\"${ocdbStorage}\")"
  if [[ -n $pretend ]]; then
    touch $qaMergedOutputFileName
    touch merge.log
    touch trending.root
  else
    #aliroot -l -b -q "merge.C(\"$qaFilesToMerge\",\"\",kFALSE,\"$qaMergedOutputFileName\")"
    aliroot -b -q "QAtrain_duo.C(\"_barrel\",${runNumber},\"$qaFilesToMerge\",1,\"${ocdbStorage}\")"
    mv QAresults_barrel.root $qaMergedOutputFileName
    mv trending_barrel.root trending.root
  fi
  
  ls -ltrh

  #copy all to output dir
  cp -p --recursive ${runpath}/* ${outputDir}
  
  #validate merge cpass1
  cd ${outputDir}
  calibDoneFile="${commonOutputPath}/merge.cpass1.run${runNumber}.done"
  touch ${calibDoneFile}
  [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" > ${calibDoneFile}
  [[ -f $qaMergedOutputFileName ]] && echo "qafile $outputDir/$qaMergedOutputFileName" >> ${calibDoneFile}
  [[ -f trending.root ]] && echo "trendingfile $outputDir/trending.root" >> ${calibDoneFile}
  echo "dir ${outputDir}" >> ${calibDoneFile}
  summarizeLogs >>  ${calibDoneFile}

  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
}

goMerge()
{
  #generic root merge using CPass1 merge.C script
  inputList=$1
  outputFile=$2  
  configFile=${3-"becnhmark.config"}
  source $configFile
  source $setupAliROOTenvInCurrentShell
  aliroot -b -q "$ALICE_ROOT/PWGPP/CalibMacros/CPass0/merge.C(\"${inputList}\",\"\",kFALSE,\"${outputFile}\")" > merge_${inputList}.log
}

goMakeSummary()
{
  configFile=$1
  outputPath=$2
  source $configFile
  
  # log filtering, script needs to take the base dir as argument
  if [[ -x $logFilteringScript ]]; then
    cd ${outputPath}
    ${logFilteringScript} $outputPath &>/dev/null
    cd -
  fi

  awk 'BEGIN {nFiles=0;nCore=0;} 
       /^calibfile/ {nFiles++;} 
       /core / {nCore++i;}
       END {print     "cpass0 produced "nFiles" calib files, "nCore" core files";}' cpass0.job*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /rec.*log OK/ {nOK++;} 
       /rec.*log BAD/ {nBAD++;} 
       END {print     "cpass0 reco:  OK: "nOK" BAD: "nBAD;}' cpass0.job*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /calib.*log OK/ {nOK++;} 
       /calib.*log BAD/ {nBAD++;} 
       END {print "cpass0 calib: OK: "nOK" BAD: "nBAD;}' cpass0.job*done
  
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /merge.log OK/ {nOK++;} 
       /merge.log BAD/ {nBAD++;} 
       END {print "cpass0 merge: OK: "nOK" BAD: "nBAD;}' merge.cpass0*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /ocdb.log OK/ {nOK++;} 
       /ocdb.log BAD/ {nBAD++;} 
       END {print   "cpass0 OCDB:  OK: "nOK" BAD: "nBAD;}' merge.cpass0*done
  
  awk 'BEGIN {nFiles=0;nCore=0;} 
       /^calibfile/ {nFiles++;} 
       /core / {nCore++;}
       END {print     "cpass1 produced "nFiles" calib files, "nCore" core files";}' cpass1.job*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /rec.*log OK/ {nOK++;} 
       /rec.*log BAD/ {nBAD++;} 
       END {print     "cpass1 reco:  OK: "nOK" BAD: "nBAD;}' cpass1.job*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /calib.*log OK/ {nOK++;} 
       /calib.*log BAD/ {nBAD++;} 
       END {print "cpass1 calib: OK: "nOK" BAD: "nBAD;}' cpass1.job*done

  awk 'BEGIN {nOK=0; nBAD=0; } 
       /merge.log OK/ {nOK++;} 
       /merge.log BAD/ {nBAD++;} 
       END {print "cpass1 merge: OK: "nOK" BAD: "nBAD;}' merge.cpass1*done
  awk 'BEGIN {nOK=0; nBAD=0; } 
       /ocdb.log OK/ {nOK++;} 
       /ocdb.log BAD/ {nBAD++;} 
       END {print   "cpass1 OCDB:  OK: "nOK" BAD: "nBAD;}' merge.cpass1*done
  
  #if set email the summary
  [[ -n $mailSummaryTo ]] && cat $log | mail -s "benchmark $productionID done" $mailSummaryTo

  return 0
}

goSubmitMakeflow()
{
  #run
  inputFileList=$1
  productionID=$2
  configFile=$3
  runNumber=$4

  [[ -z ${configFile} ]] && configFile="benchmark.config"
  [[ ! -f ${configFile} ]] && echo "no config file found (${configFile})" && return 1
  source $configFile

  if [[ ! $(which makeflow &>/dev/null) && -n $makeflowPath ]]; then
    echo "setting the makflow path from the config: "
    echo "  export PATH=${makeflowPath}:${PATH}"
    export PATH=${makeflowPath}:${PATH}
  fi

  #submit - use makeflow if available, fall back to old stuff when makeflow not there
  if which makeflow; then
    mkdir -p $productionID
    cp $0 $productionID
    cp $configFile $productionID
    cp $inputFileList $productionID
    cd $productionID
    goGenerateMakeflow "$@" > benchmark.makeflow

    [[ -n $createWorkqueue ]] && createWorkqueue=0
    if [[ $createWorkqueue -eq 1 ]]; then
      echo "creating workqueue:"
      echo "${createWorkqueueCommand}"
      eval ${createWorkqueueCommand}
    fi
    
    makeflow ${makeflowOptions} benchmark.makeflow
    cd ../
  else 
    echo "no makeflow!"
  fi
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
              "QAtrain_duo.C"
              "runCPass1.sh"
              "recCPass1.C"
              "recCPass1_OuterDet.C"
              "runCalibTrain.C"
              "runCPass0.sh"
              "recCPass0.C"
  )
  for file in ${inputFiles[*]}; do
    [[ -f ${file} ]] && copyFiles+=("${file}")
  done

  #create the makeflow file
  [[ -n ${batchFlags} ]] && echo "BATCH_OPTIONS = ${batchFlags}"
  declare -a arr_cpass1_final
  declare -a arr_cleanup
  listOfRuns=${runNumber}
  [[ -z ${runNumber} ]] && listOfRuns=($(while read x; do guessRunNumber $x; done < ${inputFileList} | sort | uniq))
  runindex=0
  for runNumber in ${listOfRuns[*]}; do
    [[ -z $runNumber ]] && continue
    [[ ! ${runNumber} =~ ^[0-9]*[0-9]$ ]] && continue

    unset arr_cpass0_outputs_onerun
    unset arr_cpass1_outputs_onerun
    declare -a arr_cpass0_outputs_onerun
    declare -a arr_cpass1_outputs_onerun

    jobindex=0
    while read inputFile; do
      currentDefaultOCDB=${defaultOCDB}
      [[ -z ${autoOCDB} ]] && autoOCDB=1
      if [[ ${autoOCDB} -ne 0 ]]; then
        currentDefaultOCDB=$(setYear $inputFile $defaultOCDB)
      fi

      #CPass0
      arr_cpass0_outputs_onerun[$jobindex]="cpass0.job${jobindex}.run${runNumber}.done"
      echo "${arr_cpass0_outputs_onerun[$jobindex]} : benchmark.sh ${configFile} ${copyFiles[@]}"
      echo " ${alirootEnv} ./benchmark.sh CPass0 ${commonOutputPath}/cpass0/000${runNumber} $inputFile $nEvents $currentDefaultOCDB $configFile $runNumber $jobindex"
      echo

      #CPass1
      arr_cpass1_outputs_onerun[$jobindex]="cpass1.job${jobindex}.run${runNumber}.done"
      echo "${arr_cpass1_outputs_onerun[$jobindex]} : benchmark.sh ${configFile} cpass0.localOCDB.${runNumber}.tgz ${copyFiles[@]}"
      echo " ${alirootEnv} ./benchmark.sh CPass1 ${commonOutputPath}/cpass1/000${runNumber} $inputFile $nEvents $currentDefaultOCDB $configFile $runNumber $jobindex"
      echo
      ((jobindex++))

    done< <(grep "/000$runNumber/" $inputFileList)
    
    #CPass0 list of Calib files to merge
    arr_cpass0_calib_list[$runindex]="cpass0.calib.run${runNumber}.list"
    echo "cpass0.calib.run${runNumber}.list: benchmark.sh ${arr_cpass0_outputs_onerun[*]}"
    echo "  ./benchmark.sh printValues calibfile cpass0.calib.run${runNumber}.list ${arr_cpass0_outputs_onerun[*]} "
    echo

    #CPass1 list of Calib/QA files to merge
    # the trick is to have the string "Stage.txt" in the file name of the list of directories with QA output to trigger
    # the production of the trending tree (only then the task->Finish() will be called in QAtrain_duo.C, on the grid
    # this corresponds to the last merging stage)
    arr_cpass1_calib_list[$runindex]="cpass1.calib.run${runNumber}.list"
    echo "cpass1.calib.run${runNumber}.list : benchmark.sh ${arr_cpass1_outputs_onerun[*]}"
    echo "  ./benchmark.sh printValues calibfile cpass1.calib.run${runNumber}.list ${arr_cpass1_outputs_onerun[*]};"
    echo
    arr_cpass1_QA_list[$runindex]="cpass1.QA.run${runNumber}.lastMergingStage.txt.list"
    echo "cpass1.QA.run${runNumber}.lastMergingStage.txt.list: benchmark.sh ${arr_cpass1_outputs_onerun[*]}"
    echo "  ./benchmark.sh printValues dir cpass1.QA.run${runNumber}.lastMergingStage.txt.list ${arr_cpass1_outputs_onerun[*]}"
    echo

    #CPass0 merging
    arr_cpass0_final[$runindex]="merge.cpass0.run${runNumber}.done"
    echo "cpass0.localOCDB.${runNumber}.tgz ${arr_cpass0_final[$runindex]}: cpass0.calib.run${runNumber}.list benchmark.sh ${configFile} ${copyFiles[@]}"
    echo " ${alirootEnv} ./benchmark.sh MergeCPass0 ${commonOutputPath}/cpass0/000${runNumber} $currentDefaultOCDB ${configFile} $runNumber cpass0.calib.run${runNumber}.list"
    echo

    #CPass1 Calib/QA merging
    arr_cpass1_final[$runindex]="merge.cpass1.run${runNumber}.done"
    echo "${arr_cpass1_final[$runindex]}: cpass0.localOCDB.${runNumber}.tgz cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.lastMergingStage.txt.list benchmark.sh ${configFile} ${copyFiles[@]}"
    echo " ${alirootEnv} ./benchmark.sh MergeCPass1 ${commonOutputPath}/cpass1/000${runNumber} $currentDefaultOCDB ${configFile} $runNumber cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.lastMergingStage.txt.list"
    echo
    
    if [[ -n $runESDfiltering ]]; then
      #CPass1 list of final ESD files
      arr_cpass1_ESD_list[$runindex]="cpass1.ESD.run${runNumber}.list"
      echo "${arr_cpass1_ESD_list[$runindex]} : benchmark.sh ${arr_cpass1_outputs_onerun[*]}"
      echo " ./benchmark.sh printValues esd ${arr_cpass1_ESD_list[$runindex]} ${arr_cpass1_outputs_onerun[*]}"
      echo

      #run the filtering
      arr_cpass1_filtering[${runindex}]="filtering.cpass1.run${runNumber}.done"
      echo "${arr_cpass1_filtering[${runindex}]} : benchmark.sh ${configFile} ${arr_cpass1_ESD_list[$runindex]} cpass0.localOCDB.${runNumber}.tgz"
      echo " LOCAL ./benchmark.sh MakeFilteredTrees ${commonOutputPath}/cpass1/000${runNumber}/ ${runNumber} ${arr_cpass1_ESD_list[$runindex]} ${filteringFactorHighPt} ${filteringFactorV0s} ${currentDefaultOCDB}"
      echo
    fi

    arr_cleanup[$runindex]="${arr_cpass0_outputs_onerun[*]} ${arr_cpass1_outputs_onerun[*]} cpass0.localOCDB.${runNumber}.tgz"

    ((runindex++))
  done

  ## final output:
  #CPass1 list of final Calib/QA files
  echo "calib.list : benchmark.sh ${arr_cpass1_final[*]}"
  echo " ./benchmark.sh printValues calibfile calib.list ${arr_cpass1_final[*]};"
  echo
  echo "qa.list : benchmark.sh ${arr_cpass1_final[*]}"
  echo " ./benchmark.sh printValues qafile qa.list ${arr_cpass1_final[*]}"
  echo
  echo "trending.list : benchmark.sh ${arr_cpass1_final[*]}"
  echo " ./benchmark.sh printValues trendingfile trending.list ${arr_cpass1_final[*]}"
  echo 
  if [[ -n $runESDfiltering ]]; then
  #list of filtered files
    echo "filtering.list : benchmark.sh ${arr_cpass1_filtering[*]}"
    echo " ./benchmark.sh printValues filteredTree filtering.list ${arr_cpass1_filtering[*]}"
    echo
  fi

  #merge trending
  echo "trending_merged.root : benchmark.sh benchmark.config trending.list"
  echo " ./benchmark.sh merge trending.list trending_merged.root ${configFile}"

  #makeQAplots
  if [[ $runMakeQAplots -eq 1 ]]; then
    echo "QAplots/ : benchmark.sh benchmark.config qa.list"
    echo " ./benchmark.sh CreateQAplots qa.list $productionID QAplots $configFile"
  fi

  #Summary
  echo "summary_${productionID}.log : ${arr_cpass0_outputs_onerun[*]} ${arr_cpass1_outputs_onerun[*]}  ${arr_cpass1_final[*]}  ${arr_cpass0_final[*]} benchmark.sh ${configFile}"
  echo " LOCAL ./benchmark.sh makeSummary ${configFile} ${commonOutputPath} |tee summary_${productionID}.log"
  echo

}

goPrintValues()
{
  #print the values given the key from any number of files (store in output file on second argument)
  key=$1
  outputFile=$2
  shift 2 #remove 2 first arguments from arg list to only pass the input files to awk
  awk -v key=${key} '$0 ~ "^"key" " {print $2}' "$@" | tee $outputFile
}

goCreateQAplots()
{
  umask 0002
  mergedQAfileList=$1
  productionID=$2
  outputDir=$3
  configFile=$4

  source $configFile
  [[ -f ${setupAliROOTenvInCurrentShell} ]] && source $setupAliROOTenvInCurrentShell

  [[ -z $logOutputDir ]] && logOutputDir=$PWD
  [[ -z $dontRedirectStdOutToLog ]] && exec 2>&1 > $logOutputDir/makeQAplots.log
  echo "$0 $*"

  [[ -z "$qaPlotsScript" ]] && echo "qaPlotsScript not defined"&&exit 1
  
  echo $qaPlotsScript "$productionID" "cpass1" $mergedQAfileList $outputDir
  $qaPlotsScript "$productionID" "cpass1" $mergedQAfileList $outputDir
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

setYear()
{
  #set the year
  #  $1 - year to be set
  #  $2 - where to set the year
  year1=$(guessYear $1)
  year2=$(guessYear $2)
  local path=$2
  [[ $year1 -ne $year2 && -n $year2 && -n $year1 ]] && path=${2/\/$year2\//\/$year1\/}
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

validateLog()
{
  log=$1
  errorConditions=(
                    "There was a crash"
                    "floating"
                    "error while loading shared libraries"
                    "std::bad_alloc"
                    "s_err_syswatch_"
                    "Thread [0-9]* (Thread"
                    "AliFatal"
                    "core dumped"
  )
  
    local logstatus=0
    local errorSummary=""
    for ((i=0; i<${#errorConditions[@]};i++)); do
      local tmp=$(grep -m1 -e "${errorConditions[$i]}" $log)
      [[ -n $tmp ]] && tmp+=" : "
      errorSummary+=$tmp
    done
    if [[ -n $errorSummary ]]; then 
      logstatus=1
      echo "$errorSummary"
    fi
  
  return $logstatus
}

summarizeLogs()
{
  #print a summary of logs
  logFiles=(
            "*.log"
            "stdout"
            "stderr"
  )

  errorConditions=(
                    "There was a crash"
                    "floating"
                    "error while loading shared libraries"
                    "std::bad_alloc"
                    "s_err_syswatch_"
                    "Thread [0-9]* (Thread"
                    "AliFatal"
                    "core dumped"
  )
  
  [[ -z ${outputDir} ]] && outputDir=$PWD
  #check logs
  logstatus=0
  for log in ${logFiles[*]}; do
    finallog=${outputDir%/}/${log}
    [[ ! -f $log ]] && continue
    errorSummary=$(validateLog $log)
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
  
  #report core files
  find ${PWD} -name core
  
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
    local localOCDBpathCPass0="${PWD}/OCDB"
  fi

  echo
  echo creating the specific storage script
  echo   localOCDBaccessConfig.C
  echo   based on OCDB: $localOCDBaccessConfig
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

goMakeFilteredTrees()
{
  outputDir=$1
  runNumber=$2
  #get path to input list
  inputListfiles=$3
  #get scale number for tracks
  filterT=$4
  #get scale number for V0s
  filterV=$5
  #get OCDB path (optional)
  OCDBpath=${6}
  #get max number of files 
  maxFiles=${7-"1000000"}
  #get offset of first file
  offsetFile=${8-"0"}
  #get max number of events
  maxEvents=${9-"30000000"}
  #get offset of first event
  offsetEvent=${10-"0"}
  configFile=${11-"benchmark.config"}
  esdFileName=${12-"AliESDs_Barrel.root"}
  source ${configFile}
  [[ -n ${setupTrustedAliROOTenvInCurrentShell} ]] && source ${setupAliROOTenvInCurrentShell}

  commonOutputPath=${PWD}
  runpath=${commonOutputPath}/rundir_filtering_run${runNumber}

  mkdir -p ${outputDir}
  mkdir -p ${runpath}

  cd ${runpath}
  
  cat > filtering.log << EOF
  goMakeFilteredTrees config:
  runpath=$runpath
  outputDir=$outputDir
  commonOutputPath=$commonOutputPath
  ALICE_ROOT=$ALICE_ROOT
  PATH=$PATH
EOF

  if [[ -z ${pretend} ]];then
    aliroot -l -b -q "$filteringRunMacro(\"$inputListfiles\",$filterT,$filterV,\"$OCDBpath\",$maxFiles,$offsetFile,$maxEvents,$offsetEvent,\"$esdFileName\")" &>> filtering.log
  else
    touch filtering.log Filtered.root
  fi
  pwd
  ls -ltrh
  echo mv -f * ${outputDir}
  mv -f * ${outputDir}
  doneFile=${commonOutputPath}/filtering.cpass1.run${runNumber}.done
  touch ${doneFile}
  [[ -f ${outputDir}/Filtered.root ]] && echo "filteredTree ${outputDir}/Filtered.root" > ${doneFile}
  [[ -f ${outputDir}/filtering.log ]] && echo "${outputDir}/filtering.log"
  cd ${commonOutputPath}
  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
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

  [[ -z $waitForJOBID ]] && waitForJOBID=0

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
    if [[ ${waitForJOBID} -eq 0 ]]; then
        echo $batchCommand $batchFlags -J "$JobID[$startID-$endID]" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
        $batchCommand $batchFlags -J "$JobID[$startID-$endID]" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
      else
        echo $batchCommand $batchFlags -J "$JobID[$startID-$endID]" -w "ended($waitForJOBID)" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
        $batchCommand $batchFlags -J "$JobID[$startID-$endID]" -w "ended($waitForJOBID)" -e "$targetDirectory/logs/job_%I.err" -o "$targetDirectory/logs/job_%I.out" "$command"     
      fi
      startID=`expr $endID + 1`
    done
  else 
    #new SGE farm
    if [[ ${waitForJOBID} =~ "000" ]]; then
      echo $batchCommand $batchFlags -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
      $batchCommand $batchFlags -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
    else
      echo $batchCommand $batchFlags -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -hold_jid "$waitForJOBID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
      $batchCommand $batchFlags -wd ${targetDirectory} -b y -V -N "$JobID" -t "$startID-$endID" -hold_jid "$waitForJOBID" -e "$targetDirectory/logs/" -o "$targetDirectory/logs/" "$command" $commandArgs
    fi
  fi
}

goSubmitBatch()
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

  #move the script, config and some other stuff to $commonOutputPath first, then use them from there
  self=$(readlink -f "$0")
  configPath=$(dirname $configFile)
  export commonOutputPath=${baseOutputDirectory}/${productionID}
  mkdir -p ${commonOutputPath}
  cp $self ${commonOutputPath}
  cp $configFile ${commonOutputPath}
  cp $inputList ${commonOutputPath}
  self=${commonOutputPath}/${self##*/}
  chmod u+x $self
  configFile=${commonOutputPath}/${configFile##*/}
  inputList=${commonOutputPath}/${inputList##*/}

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
  echo "    batchFlags:      $batchFlags"
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
  date=$(date +%Y_%m_%d_%H%M%S)
  #for each run we submit one jobarray:
  for runNumber in ${listOfRuns[*]}; do
    
    [[ -z $runNumber ]] && continue
    [[ ! ${runNumber} =~ ^[0-9]*[0-9]$ ]] && continue

    JOBpostfix="${productionID//"/"/_}_${runNumber}_${date}"
    JOBID1="p0_${JOBpostfix}"
    JOBID1wait="w0_${JOBpostfix}"
    JOBID2="m0_${JOBpostfix}"
    JOBID2wait="wm0_${JOBpostfix}"
    JOBID3="op0_${JOBpostfix}"
    JOBID3wait="wop0_${JOBpostfix}"
    JOBID4="p1_${JOBpostfix}"
    JOBID4wait="w1_${JOBpostfix}"
    JOBID5="m1_${JOBpostfix}"
    JOBID5wait="wm1_${JOBpostfix}"
    JOBID6="s1_${JOBpostfix}"
    JOBID6wait="ws1_${JOBpostfix}"
    JOBID7="QA_${JOBpostfix}"
    JOBmakeESDlistCPass1="lp1_${JOBpostfix}"
    JOBfilterESDcpass1="fp1_${JOBpostfix}"
    LASTJOB="000"

    oneInputFile=$(egrep -m1 "$runNumber\/" ${inputList})

    currentDefaultOCDB=${defaultOCDB}
    [[ -z ${autoOCDB} ]] && autoOCDB=1
    if [[ ${autoOCDB} -ne 0 ]]; then
      currentDefaultOCDB=$(setYear $oneInputFile $defaultOCDB)
    fi

    echo "submitting run $runNumber with OCDB $currentDefaultOCDB"

    ################################################################################
    ################################################################################
    # run the CPass0 if requested

    if [ $runCPass0reco -eq 1 ]; then

      echo
      echo "starting CPass0... for run $runNumber"
      echo

      # create directory and copy all files that are needed
      targetDirectory="${commonOutputPath}/000${runNumber}/cpass0"
      mkdir -p $targetDirectory
      mkdir -p $targetDirectory/logs
      filesCPass0=( 
      "$configPath/runCPass0.sh"
      "$configPath/recCPass0.C"
      "$configPath/runCalibTrain.C"
      "$configPath/localOCDBaccessConfig.C"
      "$configPath/OCDB.root"
      )
      for file in ${filesCPass0[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done


      localInputList=$targetDirectory/${inputList##*/}
      rm -f $localInputList
      egrep "\/000$runNumber\/" $inputList >> $localInputList
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

      submit $JOBID1 1 $nFiles 000 "$self" "CPass0 $targetDirectory $localInputList $nEvents $currentDefaultOCDB $configFile $runNumber"

      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "$JOBID1wait" 1 1 000 "$self" "WaitForOutput ${commonOutputPath} 'cpass0.job*.run$runNumber.done' $nFilesToWaitFor $maxSecondsToWait '-maxdepth 1'"
      LASTJOB=$JOBID1wait

    fi #end running CPass0
    ################################################################################


    ################################################################################
    # submit merging of CPass0, depends on the reconstruction

    if [ $runCPass0MergeMakeOCDB -eq 1 ]; then

      echo
      echo "submit CPass0 merging for run $runNumber"
      echo

      targetDirectory="${commonOutputPath}/000${runNumber}/cpass0"
      mkdir -p $targetDirectory

      #copy the scripts
      filesMergeCPass0=(
      "$configPath/OCDB.root"
      "$configPath/mergeMakeOCDB.byComponent.sh"
      "$configPath/mergeMakeOCDB.sh"
      "$configPath/localOCDBaccessConfig.C"
      "$configPath/mergeByComponent.C"
      "$configPath/makeOCDB.C"
      "$configPath/merge.C"
      )
      for file in ${filesMergeCPass0[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done
  

      echo submit $JOBID2 1 1 "$LASTJOB" "$self" "MergeCPass0 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      submit $JOBID2 1 1 "$LASTJOB" "$self" "MergeCPass0 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      LASTJOB=$JOBID2

      echo
    fi
    # end of merging CPass0
    ################################################################################

    ################################################################################
    ################################################################################
    # run the CPass1 if requested

    if [ $runCPass1reco -eq 1 ]; then

      targetDirectory="${commonOutputPath}/000${runNumber}/cpass1"

      # safety feature: if we are re-running for any reason we want to delete the previous output first.
      [[ -d $targetDirectory ]] && rm -rf $targetDirectory/* && echo "removed old output at $targetDirectory/*"

      echo
      echo "starting CPass1... for run $runNumber"
      echo

      # create directory and copy all files that are needed
      mkdir -p $targetDirectory
      mkdir -p $targetDirectory/logs
      
      cp -f $configFile $targetDirectory
      filesCPass1=( 
      "$configPath/runCPass1.sh"
      "$configPath/recCPass1.C"
      "$configPath/recCPass1_OuterDet.C"
      "$configPath/runCalibTrain.C"
      "$configPath/QAtrain_duo.C"
      "$configPath/localOCDBaccessConfig.C"
      "$configPath/OCDB.root"
      )
      for file in ${filesCPass1[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      localInputList=$targetDirectory/${inputList##*/}
      rm -f $localInputList
      egrep "\/000$runNumber\/" $inputList >> $localInputList
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

      submit $JOBID4 1 $nFiles "$LASTJOB" "$self" "CPass1 $targetDirectory $localInputList $nEvents $currentDefaultOCDB $configFile $runNumber"

      ################################################################################
      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "$JOBID4wait" 1 1 "$LASTJOB" "$self" "WaitForOutput ${commonOutputPath} 'cpass1.job*.run${runNumber}.done' $nFilesToWaitFor $maxSecondsToWait '-maxdepth 1'"
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

      targetDirectory="${commonOutputPath}/000${runNumber}/cpass1"
      mkdir -p $targetDirectory

      # copy files 
      filesMergeCPass1=(
      "$configPath/OCDB.root"
      "$configPath/localOCDBaccessConfig.C"
      "$configPath/mergeMakeOCDB.byComponent.sh"
      "$configPath/mergeByComponent.C"
      "$configPath/makeOCDB.C"
      "$configPath/merge.C"
      "$configPath/mergeMakeOCDB.sh"
      "$configPath/QAtrain_duo.C"
      )
      for file in ${filesMergeCPass1[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      echo submit "$JOBID5" 1 1 "$LASTJOB" "$self" "MergeCPass1 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      submit "$JOBID5" 1 1 "$LASTJOB" "$self" "MergeCPass1 $targetDirectory $currentDefaultOCDB $configFile $runNumber"
      LASTJOB=$JOBID5
      echo
    fi

    ##############################
    if [ $runESDfiltering -eq 1 ]; then
      export targetDirectory=$commonOutputPath
      echo submitting filtering
      echo targetDirectory=$targetDirectory
      submit "$JOBmakeESDlistCPass1" 1 1 "$LASTJOB" "$self" "printValues esd ${commonOutputPath}/cpass1.ESD.run${runNumber}.list ${commonOutputPath}/cpass1.job*.run${runNumber}.done "
      submit "$JOBfilterESDcpass1" 1 1 "$JOBmakeESDlistCPass1" "$self" "MakeFilteredTrees ${commonOutputPath}/000${runNumber}/cpass1 ${runNumber} ${commonOutputPath}/cpass1.ESD.run${runNumber}.list ${filteringFactorHighPt} ${filteringFactorV0s} ${currentDefaultOCDB} 1000000 0 10000000 0 ${configFile}"
    fi

  done

  #########################################
  #make lists with output files - QA, trending, filtering and calibration
  export targetDirectory=$commonOutputPath
  submit "flQA_${JOBpostfix}" 1 1 "$LASTJOB" "$self" "printValues qafile ${commonOutputPath}/qa.list ${commonOutputPath}/merge.cpass1.run*.done"
  submit "flCal_${JOBpostfix}" 1 1 "$LASTJOB" "$self" "printValues calibfile ${commonOutputPath}/calib.list ${commonOutputPath}/merge.cpass1.run*.done"
  submit "flTr_${JOBpostfix}" 1 1 "$LASTJOB" "$self" "printValues trendingfile ${commonOutputPath}/trending.list ${commonOutputPath}/merge.cpass1.run*.done"
  if [[ -n $runESDfiltering ]]; then
    submit "wFil_${JOBpostfix}" 1 1 "$LASTJOB" "$self" "WaitForOutput ${commonOutputPath} 'filtering.cpass1.run*.done' ${#listOfRuns[*]} $maxSecondsToWait '-maxdepth 1'"
    submit "flFil_${JOBpostfix}" 1 1 "wFil_${JOBpostfix}" "$self" "printValues filteredTree ${commonOutputPath}/filtering.list ${commonOutputPath}/filtering.cpass1.run*.done"
  fi

  #merge trending
  submit "mTr_${JOBpostfix}" 1 1 "flTr_${JOBpostfix}" "$self" "merge trending.list trending_merged.root ${configFile}"
  

  ################################################################################
  ################################################################################
  [[ -z $runMakeSummary ]] && runMakeSummary=0
  if [ $runMakeSummary -eq 1 ]; then
    echo
    echo "submit make a summary"
    echo

    export targetDirectory="${commonOutputPath}"
    mkdir -p $targetDirectory
    mkdir -p $targetDirectory/logs
    [[ ! -f $targetDirectory/${self##*/} ]] && cp -f $self $targetDirectory
    submit "$JOBID6" 1 1 "$LASTJOB" "$self" "makeSummary $configFile $targetDirectory $commonOutputPath"
  fi
  ################################################################################
  
  ################################################################################
  ################################################################################
  if [ $runMakeQAplots -eq 1 ]; then
    echo
    echo "submit make QA plots"
    echo

    mkdir -p ${commonOutputPath}/logs
    targetDirectory="${commonOutputPath}"

    submit "$JOBID7" 1 1 "flQA_${JOBpostfix}" "$self" "CreateQAplots qa.list $productionID QAplots $configFile"
    LASTJOB=$JOBID7
  fi
  ################################################################################

  #restore stdout
  exec 1>&7 7>&-
  echo "jobs submitted."
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
  while true; do
    #n=$(find $searchPath ${extraFindOptions} -name "$fileName" | wc -l)
    n=$(/bin/ls -1 ${searchPath}/${fileName} 2>/dev/null | wc -l)
    [[ $n -gt 0 ]] && echo "found $n X $fileName"
    [[ $n -ge $numberOfFiles ]] && break
    [[ $SECONDS -gt $maxSecondsToWait ]] && break
    sleep 60
  done
  echo "DONE! exiting..."
}

main "$@"
