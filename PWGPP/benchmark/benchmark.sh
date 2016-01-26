#!/usr/bin/env bash
#include benchmark.config

# Use extended glob. We can do negative matches for instance. This *has* to be defined here.
shopt -s extglob

# blame: Mikolaj Krzewicki, mkrzewic@cern.ch
# this script runs the CPass0/CPass1 train
# produced OCDB updates are local

#some defaults
#autoOCDB=0
defaultOCDB="raw://"
#runNumber=167123
#makeflowPath="/hera/alice/aux/cctools/bin"
#makeflowOptions="-T wq -N alice -d all -C ali-copilot.cern.ch:9097"
#makeflowOptions="-T wq -N alice -C ali-copilot.cern.ch:9097"
makeflowOptions=""
#batchCommand="/usr/bin/qsub"
batchFlags=""
baseOutputDirectory="$PWD/output"
#alirootEnv="/cvmfs/alice.cern.ch/bin/alienv setenv AliRoot/v5-04-34-AN -c"
#alirootEnv="/home/mkrzewic/alisoft/balice_master.sh"
#trustedQAtrainMacro='/hera/alice/mkrzewic/gsisvn/Calibration/QAtrain_duo.C'
reconstructInTemporaryDir=0
recoTriggerOptions="\"\""
percentProcessedFilesToContinue=100
maxSecondsToWait=$(( 3600*24 ))
nEvents=-1
nMaxChunks=0
postSetUpActionCPass0=""
postSetUpActionCPass1=""
runCPass0reco=1
runCPass0MergeMakeOCDB=1
runCPass1reco=1
runCPass1MergeMakeOCDB=1
runESDfiltering=1
filteringFactorHighPt=1e2
filteringFactorV0s=1e1
MAILTO=""
#pretend=1
#dontRedirectStdOutToLog=1
logToFinalDestination=1
ALIROOT_FORCE_COREDUMP=1
pretendDelay=0
copyInputData=0

main()
{
  #run in proper mode depending on the selection
  source $ALICE_PHYSICS/PWGPP/scripts/alilog4bash.sh false
  source $ALICE_PHYSICS/PWGPP/scripts/utilities.sh false
  if [[ $# -lt 1 ]]; then
    if [[ ! "${0}" =~ "bash" ]]; then
      echo "uses makeflow:"
      echo " ${0} \"run\" productionID inputList configFile [extraOpts]"
      echo "uses a batch system (SGE):"
      echo " ${0} \"submit\" productionID inputList configFile [extraOpts]"
      echo "extraOpts if non-empty override the config file, e.g.:"
      echo " ${0} submit test1 benchmark.list benchmark.config runNumber=169123 nEvents=10"
    fi
    return
  fi

  #define some aliases - default is to call one of the functions directly
  runMode=${1}
  umask 0002
  shift
  case ${runMode} in
    "CPass0") goCPass0 "$@";;
    "CPass1") goCPass1 "$@";;
    "MakeLocalOCDBaccessConfig") goMakeLocalOCDBaccessConfig "$@";;
    "MergeCPass0") goMergeCPass0 "$@";;
    "MergeCPass1") goMergeCPass1 "$@";;
    "MakeFilteredTrees") goMakeFilteredTrees "$@";;
    "MakeSummary") goMakeSummary "$@";;
    "run") goSubmitMakeflow "$@";;
    "submit") goSubmitBatch "$@";;
    "test") goTest "$@";;
    "GenerateMakeflow") goGenerateMakeflow "$@";;
    "PrintValues") goPrintValues "$@";;
    "CreateQAplots") goCreateQAplots "$@";;
    "WaitForOutput") goWaitForOutput "$@";;
    "Merge") goMerge "$@";;
    *) 
      ${runMode} "$@"
    ;;
  esac
  return 0
}

generateMC()
{
  #generate one raw chunk in current directory
  SEED=${JOB_ID}${SGE_TASK_ID}
  export CONFIG_SEED=${SEED}
  export runNumber=${1}
  OCDBpath=${2}
  nEventsim=${3}
  if [[ -n ${pretend} ]]; then
    sleep ${pretendDelay}
    touch galice.root
  else
    if [[ -f sim.C && -f Config.C ]] ; then
        time aliroot -b -q -x sim.C\(${runNumber},\"${OCDBpath}\",${nEventsim}\) >sim.log 2>&1
        mv syswatch.log simwatch.log
    fi
  fi
}

goCPass0() (
  # Wrapper function that calls goCPass with the CPass0 option.
  declare -a params
  params=("$1" "$2" "$3" "$4" "$5" "$6" "$7" CPass0)
  shift 7
  goCPass "${params[@]}" "$@"
)

goCPass1() (
  # Wrapper function that calls goCPass with the CPass1 option.
  declare -a params
  params=("$1" "$2" "$3" "$4" "$5" "$6" "$7" CPass1)
  shift 7
  goCPass "${params[@]}" "$@"
)

goCPass()
(
  umask 0002
  
  targetDirectory=$1
  inputList=$2
  nEvents=$3
  ocdbPath=$4
  configFile=$5
  export runNumber=$6
  jobindex=$7
  cpass=$8  # cpass=CPass0 or CPass1
  shift 8
  extraOpts=("$@")
  cpass=${cpass##*CPass}

  parseConfig configFile=$configFile "$@" || return 1

  echo Start: goCPass${cpass}
  alilog_info "[BEGIN] goCPass${cpass}() with following extra parameters $*"

  # Remember working directory provided by the batch system (i.e. current dir).
  batchWorkingDirectory=$PWD

  # Use the jobindex only if set and non-negative
  if [[ -z "$jobindex" || "$jobindex" -lt 0 ]]; then
    [[ -n "$LSB_JOBINDEX" ]] && jobindex=$LSB_JOBINDEX
    [[ -n "$SGE_TASK_ID" ]] && jobindex=$SGE_TASK_ID
    if [[ -z ${jobindex} ]]; then
      echo "no jobindex!"
      alilog_error "goCPass${cpass}() No job index [Paremeters] $*"
      return 1
    fi
  fi

  [[ -z "$commonOutputPath" ]] && commonOutputPath=$PWD

  # .done files signal job completion to Makeflow.
  [[ -n "$useProfilingCommand" ]] \
    && doneFileBase="profiling.cpass${cpass}.job${jobindex}.run${runNumber}.done" \
    || doneFileBase="cpass${cpass}.job${jobindex}.run${runNumber}.done"

  # .done file is created locally as $doneFileTmp and copied on the remote when finished.
  mkdirLocal "$commonOutputPath/meta" || return 1
  doneFileTmp="$batchWorkingDirectory/$doneFileBase"
  doneFile="$commonOutputPath/meta/$doneFileBase"

  [[ -f "$alirootSource" && -z "$ALICE_ROOT" ]] && source ${alirootSource}
  
  [[ -n "$ALIROOT_FORCE_COREDUMP" ]] && export ALIROOT_FORCE_COREDUMP && ulimit -c unlimited

  # The contents of this is stored in the tree and used later (e.g. AliAnalysisTaskPIDResponse)!
  # At the QA stage the pass number is guessed from the path stored here.
  # The Format is:
  #   Packages= ;OutputDir= ;LPMPass= ;TriggerAlias= ;LPMRunNumber= ;LPMProductionType= ;
  #   LPMInteractionType= ;LPMProductionTag= ;LPMAnchorRun= ;LPMAnchorProduction= ;LPMAnchorYear= 
  export PRODUCTION_METADATA="OutputDir=cpass${cpass}"

  # Check if input list exists. Works both for local and remote cases.
  if ! statRemote "$inputList" && [[ -z $pretend ]]; then
    touch $doneFileTmp
    echo "Cannot find input list $inputList, exiting..." >> $doneFileTmp
    copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
    return 1
  fi

  [[ "$inputList" =~ \.root$ ]] && infile=$inputList \
                                || infile=$(sed -ne "${jobindex}p" $inputList | egrep '\s*\w*/\w*')
  chunkName=${infile##*/}

  outputDir="${targetDirectory}/${jobindex}_${chunkName%.*}"
  case "$reconstructInTemporaryDir" in
    1) runpath=$(mktemp -d -t cpass${cpass}.XXXXXX) ;;
    2) runpath=${PWD}/rundir_cpass${cpass}_${runNumber}_${jobindex} ;;
    *) runpath=$outputDir ;;
  esac

  # TODO: check if this is really for CPass1 only. Asymmetry found when merging CPass0/1.
  # TODO: Used only in MC. Check if this still works.
  if [[ $cpass == 1 ]]; then
    if [[ "${infile}" =~ galice\.root ]]; then
      ln -s ${inputList%/*}/* ${runpath}
      infile=""
    fi
  fi

  # Robust error check in directory creation. After this block we are in $runpath.
  if ! mkdirLocal "$outputDir" || ! mkdirLocal "$runpath" || ! cd "$runpath"; then
    touch $doneFileTmp
    echo "Error creating $outputDir or runpath $runpath, or cd'ing to it" >> $doneFileTmp
    copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
    return 1
  fi

  # runCPassX/C expects the raw chunk to be linked in the run dir despite it being accessed by the
  # full path.
  [[ $copyInputData == 0 ]] && ln -s ${infile} ${runpath}/${chunkName} \
                            || copyFileFromRemote ${infile} ${runpath}/

  logOutputDir=$runpath
  [[ -n "$logToFinalDestination" ]] && logOutputDir=${outputDir}
  [[ -z "$dontRedirectStdOutToLog" ]] && exec &> ${logOutputDir}/stdout
  echo "$0 $*"

  ##### MC -- TODO: does this still work? is it really needed during CPass0 only?
  if [[ $cpass == 0 && -n $generateMC ]]; then
    olddir=${PWD}
    outputDirMC=${commonOutputPath}/000${runNumber}/sim/${jobindex}
    simrunpath=${outputDirMC}
    #[[ ${simulateInTemporaryDir} -eq 1 && -n ${TMPDIR} ]] && simrunpath=${TMPDIR}
    [[ ${simulateInTemporaryDir} -eq 1 ]] && simrunpath=$(mktemp -d -t cpass0MC.XXXXXX)
    mkdir -p ${outputDirMC}
    mkdir -p ${simrunpath}
    if cd ${simrunpath}; then

      filesMC=( 
      "${batchWorkingDirectory}/sim.C"
      "${batchWorkingDirectory}/rec.C"
      "${batchWorkingDirectory}/Config.C"
      "${batchWorkingDirectory}/OCDB_*.root"
      )
      for file in ${filesMC[*]}; do
        [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
      done

      generateMC ${runNumber} ${ocdbPath} ${nEvents}

      [[ ! "${simrunpath}" =~ "${outputDirMC}" ]] && mv * ${outputDirMC} #TODO check if it works
      cd ${olddir}

      ln -s ${outputDirMC}/* ${runpath}/ 

      inputList=${outputDirMC}/galice.root #TODO not valid outside shell !!!
      infile=""
    fi
  fi
  ###### /MC

  echo "#####################"
  echo CPass${cpass}:
  echo JOB setup
  echo nEvents               ${nEvents}
  echo runNumber             ${runNumber}
  echo ocdbPath              ${ocdbPath}
  echo infile                ${infile}
  echo chunkName             ${chunkName}
  echo jobindex              ${jobindex}
  echo recoTriggerOptions    ${recoTriggerOptions}
  echo targetDirectory       ${targetDirectory}
  echo commonOutputPath      ${commonOutputPath}
  echo doneFile              ${doneFile}
  echo batchWorkingDirectory ${batchWorkingDirectory}
  echo runpath               ${runpath}  
  echo outputDir             ${outputDir}
  echo PWD                   ${PWD}
  echo ALICE_ROOT            ${ALICE_ROOT}
  echo ALICE_PHYSICS         ${ALICE_PHYSICS}
  echo "########## ###########"

  alirootInfo > ALICE_ROOT.log

  # Those files are locally found. They should be copied to the current working
  # directory (which is $runpath).
  case $cpass in
    0) filesCPass=( "${batchWorkingDirectory}/runCPass0.sh"
                    "${batchWorkingDirectory}/recCPass0.C"
                    "${batchWorkingDirectory}/runCalibTrain.C"
                    "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                    "${batchWorkingDirectory}/OCDB.root"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/runCPass0.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/recCPass0.C" 
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/runCalibTrain.C" ) ;;

    1) filesCPass=( "${batchWorkingDirectory}/runCPass1.sh"
                    "${batchWorkingDirectory}/recCPass1.C"
                    "${batchWorkingDirectory}/recCPass1_OuterDet.C"
                    "${batchWorkingDirectory}/runCalibTrain.C"
                    "${batchWorkingDirectory}/QAtrain_duo.C"
                    "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                    "${batchWorkingDirectory}/${configFile}"
                    "${commonOutputPath}/meta/cpass0.localOCDB.${runNumber}.tgz"
                    "${batchWorkingDirectory}/OCDB.root"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/runCPass1.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/recCPass1.C" ) ;;
  esac

  for file in ${filesCPass[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && printExec cp -f $file . \
                                         && [[ ${file##*/} =~ .*\.sh ]] \
                                         && printExec chmod +x ${file##*/}
  done

  echo "Contents of current directory ($PWD) before running CPass${cpass}:" ; /bin/ls ; echo

  # Monkey patching: emove spaces from around arguments to root macros. For example this sometimes
  # is known to fail: root 'macro.C(argument1, argument2)'
  sed -i '/.*root .*\.C/ s|\s*,\s*|,|g' *.sh

  # A post-setup action (different for CPass0/1) can be defined. Here it is executed.
  case $cpass in
    0) postSetUpActionCPass="$postSetUpActionCPass0" ;;
    1) postSetUpActionCPass="$postSetUpActionCPass1" ;;
  esac
  [[ -n "$postSetUpActionCPass" ]] && printExec eval $postSetUpActionCPass

  # Run CPass0 or CPass1.
  case $cpass in

    0)
      # Run CPass0. Since $infile is local, by convention it must start with a slash. The initial slash
      # is there just to signal runCPass0.sh that the file is local, and it is not an actual path.
      if [[ -n "$pretend" ]]; then
        echo "Pretending to run: $runpath/runCPass0.sh /$infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
        sleep $pretendDelay
        for fakeOutput in AliESDs.root AliESDfriends.root AliESDfriends_v1.root rec.log calib.log; do
          touch $fakeOutput
        done
      else
        # Caveat: in the local case, first arg must start with a slash.
        printExec ./runCPass0.sh /$infile $nEvents $runNumber $ocdbPath $recoTriggerOptions
      fi

      # End of CPass0.
    ;;

    1)
      # This is CPass1.

      # Procure OCDB.
      echo "Downloading OCDB produced during CPass0 for run $runNumber"
      copyFileFromRemote $commonOutputPath/meta/cpass0.localOCDB.${runNumber}.tgz $PWD

      # If OCDB is found here, then create a macro that configures local OCDB access.
      # This step also decompresses the tarball into $PWD/OCDB.
      [[ -f cpass0.localOCDB.${runNumber}.tgz ]] \
        && printExec goMakeLocalOCDBaccessConfig cpass0.localOCDB.${runNumber}.tgz \
        || echo "WARNING: file cpass0.localOCDB.${runNumber}.tgz not found!"

      # Check if CPass0 has produced calibration.
      if [[ ! $(/bin/ls -1 OCDB/*/*/*/*.root 2>/dev/null) ]]; then
        touch $doneFileTmp
        echo "CPass0 produced no calibration! Exiting..." >> $doneFileTmp
        copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
        return 1
      fi

      # Create the Barrel and OuterDet directories for CPass1 and link the local OCDB directory
      # there to make the localOCDBaccessConfig.C file work, since it may point to the OCDB
      # entries using a relative path, e.g. local://./OCDB.
      echo "Linking the OCDB/ for Barrel and OuterDet directories"
      mkdir Barrel OuterDet
      ls -l
      ln -s ../OCDB Barrel/OCDB
      ln -s ../OCDB OuterDet/OCDB

      # Setup the filtering.
      # The following option enables the filtering task inside the QAtrain_duo.C
      [[ -n $runESDfiltering ]] && export QA_TaskFilteredTree=1
      #set the downscaling factors during the filtering fro expert QA (overrides the previous values)
      [[ -n ${filteringFactorHighPt} ]] && export AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF=${filteringFactorHighPt}
      [[ -n ${filteringFactorV0s} ]] && export AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF=${filteringFactorV0s}

      # Run CPass1.
      if [[ -n "$pretend" ]]; then
        echo "Pretending to run: ./runCPass1.sh /$infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
        sleep $pretendDelay
        for fakeOutput in AliESDs_Barrel.root AliESDfriends_Barrel.root AliESDfriends_v1.root \
                          QAresults_barrel.root EventStat_temp_barrel.root AODtpITS.root \
                          AliESDs_Outer.root AliESDfriends_Outer.root QAresults_outer.root \
                          EventStat_temp_outer.root rec.log calib.log qa.log filtering.log \
                          FilterEvents_Trees.root ; do
          touch $fakeOutput
        done
      else
        # Caveat: in the local case, first arg must start with a slash.
        printExec ./runCPass1.sh /$infile $nEvents $runNumber $ocdbPath $recoTriggerOptions
      fi

      # Fix output locations.
      [[ ! -f AliESDs_Barrel.root && -f Barrel/AliESDs.root ]] && mv Barrel/AliESDs.root AliESDs_Barrel.root
      [[ ! -f AliESDfriends_Barrel.root && -f Barrel/AliESDfriends.root ]] && mv Barrel/AliESDfriends.root AliESDfriends_Barrel.root
      [[ ! -f AliESDfriends_v1.root && -f Barrel/AliESDfriends_v1.root ]] && mv Barrel/AliESDfriends_v1.root .
      [[ ! -f QAresults_barrel.root && -f Barrel/QAresults_barrel.root ]] && mv Barrel/QAresults_barrel.root .
      [[ ! -f AliESDs_Outer.root && -f OuterDet/AliESDs.root ]] && mv OuterDet/AliESDs.root AliESDs_Outer.root
      [[ ! -f AliESDfriends_Outer.root && -f OuterDet/AliESDfriends.root ]] && mv OuterDet/AliESDfriends.root AliESDfriends_Outer.root
      [[ ! -f QAresults_outer.root && -f OuterDet/QAresults_outer.root ]] && mv OuterDet/QAresults_outer.root .
      [[ ! -f FilterEvents_Trees.root && -f Barrel/FilterEvents_Trees.root ]] && mv Barrel/FilterEvents_Trees.root .

      # Make the filtered tree (if requested and not already produced by QA.
      # TODO: check if it works with non-shared filesystems.
      # TODO: check if they are still needed. In principle we could remove them.
      [[ -f AliESDs_Barrel.root ]] && echo AliESDs_Barrel.root > filtered.list
      if [[ -n $runESDfiltering && ! -f FilterEvents_Trees.root && -f filtered.list ]]; then
        goMakeFilteredTrees $PWD $runNumber "$PWD/filtered.list" $filteringFactorHighPt \
                            $filteringFactorV0s $ocdbPath 1000000 0 10000000 0 \
                            $configFile AliESDs_Barrel.root "${extraOpts[@]}" > filtering.log
      fi

      # End of CPass1.
    ;;

  esac

  echo "Contents of current directory ($PWD) after running CPass$cpass:" ; /bin/ls ; echo

  # CPass has completed. Copy all created files to the destination (which might be remote).
  # Note: stdout is not copied, this will happen at the very end.
  printExec rm -f ./$chunkName
  while read cpdir; do
    printExec copyFileToRemote $cpdir/!(stdout) $outputDir/$cpdir
  done < <(find . -type d)

  # Validate CPass.
  case $cpass in

    0)
      # Validate CPass0.
      if summarizeLogs | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
        [[ -f $outputDirMC/galice.root ]] && echo "sim $outputDirMC/galice.root" >> $doneFileTmp
        [[ -f AliESDfriends_v1.root ]] && echo "calibfile $outputDir/AliESDfriends_v1.root" >> $doneFileTmp
        [[ -f AliESDs.root ]] && echo "esd $outputDir/AliESDs.root" >> $doneFileTmp
      fi
      # End of CPass0 validation.
    ;;

    1)
      # Validate CPass1.
      if summarizeLogs | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
        [[ -f AliESDs_Barrel.root ]] && echo "esd $outputDir/AliESDs_Barrel.root" >> $doneFileTmp
        [[ -f AliESDfriends_v1.root ]] && echo "calibfile $outputDir/AliESDfriends_v1.root" >> $doneFileTmp
        [[ -f QAresults_Barrel.root ]] && echo "qafile $outputDir/QAresults_Barrel.root" >> $doneFileTmp
        [[ -f QAresults_Outer.root ]] && echo "qafile $outputDir/QAresults_Outer.root" >> $doneFileTmp
        [[ -f QAresults_barrel.root ]] && echo "qafile $outputDir/QAresults_barrel.root" >> $doneFileTmp
        [[ -f QAresults_outer.root ]] && echo "qafile $outputDir/QAresults_outer.root" >> $doneFileTmp
        [[ -f FilterEvents_Trees.root ]] && echo "filteredTree $outputDir/FilterEvents_Trees.root" >> $doneFileTmp
      else
        if grep -q qa_outer.log.*OK $doneFileTmp; then
          [[ -f QAresults_Outer.root ]] && echo "qafile $outputDir/QAresults_Outer.root" >> $doneFileTmp
          [[ -f QAresults_outer.root ]] && echo "qafile $outputDir/QAresults_outer.root" >> $doneFileTmp
        fi
        if grep -q qa_barrel.log.*OK $doneFileTmp; then
          [[ -f QAresults_Barrel.root ]] && echo "qafile $outputDir/QAresults_Barrel.root" >> $doneFileTmp
          [[ -f QAresults_barrel.root ]] && echo "qafile $outputDir/QAresults_barrel.root" >> $doneFileTmp
        fi
        if grep -q filtering.log.*OK $doneFileTmp; then
          [[ -f FilterEvents_Trees.root ]] && echo "filteredTree $outputDir/FilterEvents_Trees.root" >> $doneFileTmp
        fi
      fi
      # End of CPass1 validation.
    ;;

  esac
  echo "dir $outputDir" >> $doneFileTmp

  # Copy stdout to destination.
  [[ -z "$dontRedirectStdOutToLog" ]] && echo "Copying stdout. NOTE: this is the last bit you will see in the log!"
  copyFileToRemote stdout $outputDir

  # Final cleanup (only if we are not writing directly to destination).
  [[ "$runpath" != "$outputDir" ]] && printExec rm -rf $runpath
  copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
  echo End: goCPass${cpass}
  alilog_info "[END] goCPass${cpass}() with following parameters $*"
  return 0
)

goMergeCPass0() (
  goMergeCPass CPass0 "$@"
)

goMergeCPass1() (
  goMergeCPass CPass1 "$@"
)

goMergeCPass()
(
  # Find CPass0 output files and merge them. OCDB is created.

  # First argument can be either CPass0 or CPass1.
  cpass=$1
  cpass=${cpass##*CPass}
  [[ $cpass != 0 && $cpass != 1 ]] \
    && echo "FATAL: call goMergeCPass with CPass[0|1] as first param!" && return 1
  shift

  # Arguments in common between MergeCPass0 and MergeCPass1.
  outputDir=$1
  ocdbStorage=$2
  configFile=$3
  export runNumber=$4
  calibrationFilesToMerge=$5  #can be a non-existent file, will then be produced on the fly
  shift 5

  # MergeCPass1 takes two arguments more.
  if [[ $cpass == 1 ]]; then
    qaFilesToMerge=$1
    filteredFilesToMerge=$2
    shift 2
  fi

  parseConfig configFile=$configFile "$@" || return 1

  echo Start: goMergeCPass${cpass}
  alilog_info "[BEGIN] goMergeCPass${cpass}() with following parameters $*"

  # Record the working directory provided by the batch system.
  batchWorkingDirectory=$PWD

  [[ -z "$commonOutputPath" ]] && commonOutputPath=$PWD

  # This file signals that everything went fine
  doneFileBase="merge.cpass${cpass}.run${runNumber}.done"

  # We will have two copies of the file. The tmp one is used locally, then it is
  # transferred remotely as the other one.
  doneFileTmp="${batchWorkingDirectory}/${doneFileBase}"
  doneFile="${commonOutputPath}/meta/${doneFileBase}"

  umask 0002
  ulimit -c unlimited 

  [[ -f "$alirootSource" && -z "$ALICE_ROOT" ]] && source "$alirootSource"

  case "$reconstructInTemporaryDir" in
    1) runpath=$(mktemp -d -t mergeCPass${cpass}.XXXXXX) ;;
    2) runpath=${PWD}/rundir_mergeCPass${cpass}_${runNumber} ;;
    *) runpath=$outputDir ;;
  esac

  # Robust error check in directory creation. After this block we are in $runpath.
  if ! mkdirLocal "$runpath" || ! cd "$runpath"; then
    touch "$doneFileTmp"
    echo "Error creating runpath $runpath, or cd'ing to it" >> $doneFileTmp
    copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")"
    return 1
  fi

  logOutputDir=${runpath}
  [[ -n "$logToFinalDestination" ]] && logOutputDir=${outputDir}
  [[ -z "$dontRedirectStdOutToLog" ]] && exec &> ${logOutputDir}/stdout
  echo "$0 $*"

  mergingScript="mergeMakeOCDB.byComponent.sh"

  if [[ $cpass == 1 ]]; then
    calibrationOutputFileName='AliESDfriends_v1.root'
    qaOutputFileName='QAresults*.root'
    # Important to have the string "Stage.txt" in the filename to trigger the merging.
    # It has to be a list of directories containing the files.
    qaMergedOutputFileName="QAresults_merged.root"
  fi

  echo goMergeCPass${cpass} SETUP:
  echo runNumber=${runNumber}
  echo outputDir=${outputDir}
  echo ocdbStorage=${ocdbStorage}
  echo calibrationFilesToMerge=${calibrationFilesToMerge}
  echo mergingScript=${mergingScript}
  echo commonOutputPath=${commonOutputPath}
  echo runpath=${runpath}
  [[ $cpass == 1 ]] && echo qaFilesToMerge=$qaFilesToMerge
  [[ $cpass == 1 ]] && echo calibrationOutputFileName=$calibrationOutputFileName
  
  # Copy files in case they are not already there.
  case $cpass in
    0) filesMergeCPass=( "${batchWorkingDirectory}/${calibrationFilesToMerge}"
                         "${batchWorkingDirectory}/OCDB.root"
                         "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/mergeMakeOCDB.byComponent.sh"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/mergeByComponent.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/makeOCDB.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/merge.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/mergeMakeOCDB.sh" ) ;;

    1) filesMergeCPass=( "${batchWorkingDirectory}/${calibrationFilesToMerge}"
                         "${batchWorkingDirectory}/${qaFilesToMerge}"
                         "${batchWorkingDirectory}/OCDB.root"
                         "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                         "${commonOutputPath}/meta/cpass0.localOCDB.${runNumber}.tgz"
                         "${batchWorkingDirectory}/QAtrain_duo.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.byComponent.sh"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeByComponent.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/makeOCDB.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/merge.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeMakeOCDB.sh"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/QAtrain_duo.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeQAgroups.C"
                         "${trustedQAtrainMacro}" ) ;;
  esac

  for file in ${filesMergeCPass[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && printExec cp -f $file . \
                                         && [[ ${file##*/} =~ .*\.sh ]] \
                                         && printExec chmod +x ${file##*/}
  done

  # Monkey patching: remove spaces from around arguments to root macros. For example this sometimes
  # is known to fail: root 'macro.C(argument1, argument2)'
  sed -i '/.*root .*\.C/ s|\s*,\s*|,|g' *.sh

  alirootInfo > ALICE_ROOT.log

  # Configure local OCDB storage (a macro is produced). Only used in MergeCPass1 but could be used
  # in CPass0 too. It's harmless in any case.
  if [[ -f cpass0.localOCDB.${runNumber}.tgz ]]; then
    printExec goMakeLocalOCDBaccessConfig "cpass0.localOCDB.${runNumber}.tgz"
  else
    echo "WARNING: file cpass0.localOCDB.${runNumber}.tgz not found!"
  fi

  echo "Contents of current directory before running MergeCPass${cpass} ($PWD):" ; /bin/ls ; echo
  echo "Contents of parent directory before running MergeCPass${cpass} ($(cd ..;pwd)):" ; /bin/ls .. ; echo

  # Merge calibration.
  chmod u+x $mergingScript
  mkdir -p OCDB

  # Run Merge of this CPass.
  case $cpass in
    0)
      # MergeCPass0.

      if [[ -n ${pretend} ]]; then
        sleep $pretendDelay
        for file in CalibObjects.root ocdb.log merge.log dcsTime.root; do
          touch $file
        done
        mkdir -p OCDB/TPC/Calib/{TimeGain,TimeDrift}
        echo "some calibration" >> OCDB/TPC/Calib/TimeGain/someCalibObject_0-999999_cpass0.root
        echo "some calibration" >> OCDB/TPC/Calib/TimeDrift/otherCalibObject_0-999999_cpass0.root
      else
        printExec ./$mergingScript $calibrationFilesToMerge \
                                   $runNumber \
                                   "local://./OCDB" \
                                   defaultOCDB=$ocdbStorage \
                                   fileAccessMethod=nocopy >> mergeMakeOCDB.log
    
        #produce the calib trees for expert QA (dcsTime.root)
        goMakeLocalOCDBaccessConfig ./OCDB
        printExec aliroot -b -q "${ALICE_PHYSICS}/PWGPP/TPC/macros/CalibSummary.C($runNumber,\"$ocdbStorage\")"
      fi

      # End of MergeCPass0.
    ;;

    1)
      # MergeCPass1.

      if [[ -n ${pretend} ]]; then
        sleep ${pretendDelay}
        for file in ocdb.log ${qaMergedOutputFileName} \
                    merge.log trending.root FilterEvents_Trees.root CalibObjects.root \
                    dcsTime.root; do
          touch $file
        done
        mkdir -p OCDB/TPC/Calib/{TimeGain,TimeDrift}
        echo "some calibration" >> OCDB/TPC/Calib/TimeGain/someCalibObject_0-999999_cpass1.root
        echo "some calibration" >> OCDB/TPC/Calib/TimeDrift/otherCalibObject_0-999999_cpass1.root
      else
        printExec ./${mergingScript} ${calibrationFilesToMerge} \
                                     ${runNumber} \
                                     "local://./OCDB" \
                                     defaultOCDB=${ocdbStorage} \
                                     fileAccessMethod=nocopy

        # Merge QA (and filtered trees).

        [[ -n ${AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF
        [[ -n ${AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF} ]] && export AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF

        # TODO can we delete the following line?
        #printExec aliroot -l -b -q "merge.C(\"${qaFilesToMerge}\",\"\",kFALSE,\"${qaMergedOutputFileName}\")"
        printExec aliroot -b -q "QAtrain_duo.C(\"_barrel\",${runNumber},\"${qaFilesToMerge}\",1,\"${ocdbStorage}\")" > mergeQA.log

        mv QAresults_barrel.root ${qaMergedOutputFileName}
        mv trending_barrel.root trending.root
 
        # Merge filtered trees.
        printExec aliroot -l -b -q "merge.C(\"${filteredFilesToMerge}\",\"\",kFALSE,\"FilterEvents_Trees.root\")" > mergeFilteredTrees.log

        # Produce the calib trees for expert QA.
        printExec aliroot -b -q "${ALICE_PHYSICS}/PWGPP/TPC/macros/CalibSummary.C(${runNumber},\"${ocdbStorage}\")" > calibTree.log

      fi

      # End of MergeCPass1.
    ;;
  esac
  
  # Create tarball with OCDB, store on the shared directory, create signal file on batch directory
  baseTar="cpass${cpass}.localOCDB.${runNumber}.tgz"
  printExec tar czvvf ${batchWorkingDirectory}/${baseTar} ./OCDB && \
    echo "ocdbTarball $commonOutputPath/meta/$baseTar" >> $doneFileTmp

  echo "Contents (recursive) of batch working directory after tarball creation in MergeCPass0 ($batchWorkingDirectory):" ; /bin/ls -R $batchWorkingDirectory ; echo

  # Copy all to output dir. Preserve dir structure. stdout will be copied later.
  while read cpdir; do
    copyFileToRemote $cpdir/!(stdout) $outputDir/$cpdir
  done < <(find . -type d)

  # Copy OCDB to meta.
  copyFileToRemote ${batchWorkingDirectory}/${baseTar} $commonOutputPath/meta

  # TODO: check if MC still works. Not used at CERN release validation.
  if [[ $cpass == 0 && -n ${generateMC} ]]; then
    goPrintValues sim ${commonOutputPath}/meta/sim.run${runNumber}.list ${commonOutputPath}/meta/cpass0.job*.run${runNumber}.done
  fi

  if summarizeLogs | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
    [[ -f CalibObjects.root ]] && echo "calibfile ${outputDir}/CalibObjects.root" >> ${doneFileTmp}
    [[ -f dcsTime.root ]] && echo "dcsTree ${outputDir}/dcsTime.root" >> ${doneFileTmp}

    # Those files should be there only when merging CPass1.
    [[ -f ${qaMergedOutputFileName} ]] && echo "qafile ${outputDir}/${qaMergedOutputFileName}" >> ${doneFileTmp}
    [[ -f trending.root ]] && echo "trendingfile ${outputDir}/trending.root" >> ${doneFileTmp}
    [[ -f FilterEvents_Trees.root ]] && echo "filteredTree ${outputDir}/FilterEvents_Trees.root" >> ${doneFileTmp}
  else
    grep -q "mergeQA.log.*OK" $doneFileTmp \
      && [[ -f $qaMergedOutputFileName ]] \
      && echo "qafile ${outputDir}/${qaMergedOutputFileName}" >> $doneFileTmp
    grep -q "mergeFilteredTrees.log.*OK" $doneFileTmp \
      && [[ -f FilterEvents_Trees.root ]] \
      && echo "filteredTree ${outputDir}/FilterEvents_Trees.root" >> $doneFileTmp
  fi
  echo "dir $outputDir" >> $doneFileTmp

  # Copy stdout to destination.
  [[ -z "$dontRedirectStdOutToLog" ]] && echo "Copying stdout. NOTE: this is the last bit you will see in the log!"
  copyFileToRemote stdout $outputDir

  # Final cleanup.
  [[ "$runpath" != "$outputDir" ]] && rm -rf ${runpath}
  copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
  echo End: goMergeCPass${cpass}
  alilog_info "[END] goMergeCPass${cpass}() with following parameters $*"
  return 0
)

goMerge()
(
  #generic root merge using CPass1 merge.C script
  inputList=${1}
  outputFile=${2}  
  configFile=${3-"benchmark.config"}
  shift 3
  alilog_info  "[BEGIN] goMerge() with following parameters $*"
  if ! parseConfig configFile=${configFile} "$@"; then return 1; fi
  
  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  [[ ! -f ${inputList} ]] && echo "inputList ${inputList} does not exist!" && return 1
  [[ ! -f ${configFile} ]] && echo "configFile ${configFile} does not exist!" && return 1
  umask 0002
  [[ -f ${alirootSource} && -z ${ALICE_ROOT} ]] && source ${alirootSource}
  rm -f ${outputFile}
  aliroot -b -q "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/merge.C(\"${inputList}\",\"\",kFALSE,\"${outputFile}\")" > merge_${inputList}.log
  alilog_info  "[END] goMerge() with following parameters $*"
  return 0
)

goSubmitMakeflow()
{
  #run
  productionID=${1}
  inputFileList=${2}
  configFile=${3}
  shift 3
  extraOpts=("$@")
  if ! parseConfig configFile=${configFile} "${extraOpts[@]}"; then return 1; fi
 
  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  [[ -z ${configFile} ]] && configFile="benchmark.config"
  [[ ! -f ${configFile} ]] && echo "no config file found (${configFile})" && return 1

  if [[ ! $(which makeflow &>/dev/null) && -n ${makeflowPath} ]]; then
    echo "setting the makflow path from the config: "
    echo "  export PATH=${makeflowPath}:${PATH}"
    export PATH=${makeflowPath}:${PATH}
  fi

  # Create the common output dir and the meta dir, but only if local.
  commonOutputPath=${baseOutputDirectory}/${productionID}
  mkdirLocal $commonOutputPath/meta
 
  self=$0
  
  # For reference copy the setup to the output dir.
  copyFileToRemote $self $commonOutputPath
  copyFileToRemote $configFile $commonOutputPath
  copyFileToRemote $inputFileList $commonOutputPath

  # Submit - use makeflow if available, fall back to old stuff when makeflow not there.
  if which makeflow; then
    goGenerateMakeflow $productionID $inputFileList \
                       $configFile "${extraOpts[@]}" \
                       commonOutputPath=$commonOutputPath > benchmark.makeflow
    copyFileToRemote benchmark.makeflow $commonOutputPath
    makeflow $makeflowOptions benchmark.makeflow
  else 
    echo "no makeflow!"
  fi
  
  # Summarize the run based on the makeflow log, append it to summary.log.
  awk '/STARTED/   {startTime=$3} 
       /COMPLETED/ {endTime=$3} 
       END         {print "makeflow running time: "(endTime-startTime)/1000000/3600" hours"}' \
      benchmark.makeflow.makeflowlog | tee -a summary.log
  copyFileToRemote summary.log $commonOutputPath

  return 0
}

goGenerateMakeflow()
(
  #generate the makeflow file
  [[ $# -lt 3 ]] && echo "args: id inputFileList configFile" && return 1
  productionID=${1}
  inputFileList=${2}
  configFile=${3}
  shift 3
  extraOpts=("$@")
  
  #batch systems/makeflow sometimes handle spaces in arguments poorly, so encode them
  for (( i=0;i<${#extraOpts[@]};i++ )); do 
    extraOpts[i]=$(encSpaces "${extraOpts[i]}")
  done
  extraOpts+=("encodedSpaces=1")

  if ! parseConfig configFile=${configFile} "${extraOpts[@]}" &>/dev/null; then return 1; fi
 
  #extra safety
  if [[ -z ${commonOutputPath} ]]; then
    commonOutputPath=${baseOutputDirectory}/${productionID}
    extraOpts=( "${extraOpts[@]}" "commonOutputPath=${commonOutputPath}" )
  fi

  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  [[ -z ${configFile} ]] && configFile="benchmark.config"
  [[ ! -f ${configFile} ]] && echo "no config file found (${configFile})" && return 1

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
              "runQA.sh"
	      "mergeQAgroups.C"
  )
  for file in ${inputFiles[*]}; do
    [[ -f ${file} ]] && copyFiles+=("${file}")
  done

  #create the makeflow file
  [[ -n ${batchFlags} ]] && echo "BATCH_OPTIONS = ${batchFlags}"
  declare -A arr_cpass0_merged arr_cpass1_merged
  declare -A arr_cpass0_calib_list arr_cpass1_calib_list 
  declare -A arr_cpass1_QA_list arr_cpass1_ESD_list arr_cpass1_filtered_list
  declare -A arr_cpass0_profiled_outputs
  declare -A listOfRuns
  [[ -n ${runNumber} ]] && listOfRuns[${runNumber}]=1
  while read x; do tmpRun=$(guessRunNumber ${x}); [[ -n ${tmpRun} ]] && listOfRuns[${tmpRun}]=1; done < ${inputFileList}
  for runNumber in "${!listOfRuns[@]}"; do
    [[ -z ${runNumber} ]] && continue
    [[ ! ${runNumber} =~ ^[0-9]*[0-9]$ ]] && continue

    unset arr_cpass0_outputs
    unset arr_cpass1_outputs
    declare -a arr_cpass0_outputs
    declare -a arr_cpass1_outputs

    #Header
    echo "### Automatically generated on $(LANG=C date) ###"
    echo ; echo

    jobindex=0
    inputFile=""
    while read inputFile; do
      currentDefaultOCDB=${defaultOCDB}
      [[ -z ${autoOCDB} ]] && autoOCDB=1
      if [[ ${autoOCDB} -ne 0 ]]; then
        currentDefaultOCDB=$(setYear ${inputFile} ${defaultOCDB})
      fi
      guessRunData ${inputFile}

      #Set variables
      echo "### Variables ###"
      echo "OUTPATH=\"${commonOutputPath}/${year}/${period}\""
      echo ; echo

      #CPass0
      #arr_cpass0_outputs[${jobindex}]="${commonOutputPath}/meta/cpass0.job${jobindex}.run${runNumber}.done"
      arr_cpass0_outputs[${jobindex}]="cpass0.job${jobindex}.run${runNumber}.done"
      echo "### CPass0 ###"
      echo "${arr_cpass0_outputs[${jobindex}]}: benchmark.sh ${configFile} ${copyFiles[@]}"
      echo -e "\t${alirootEnv} ./benchmark.sh CPass0 \$OUTPATH/000${runNumber}/cpass0 ${inputFile} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]}"" "
      echo ; echo

      #CPass1
      #arr_cpass1_outputs[${jobindex}]="${commonOutputPath}/meta/cpass1.job${jobindex}.run${runNumber}.done"
      arr_cpass1_outputs[${jobindex}]="cpass1.job${jobindex}.run${runNumber}.done"
      echo "### CPass1 ###"
      echo "${arr_cpass1_outputs[${jobindex}]}: benchmark.sh ${configFile} merge.cpass0.run${runNumber}.done ${copyFiles[@]}"
      echo -e "\t${alirootEnv} ./benchmark.sh CPass1 \$OUTPATH/000${runNumber}/cpass1 ${inputFile} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]}"" "
      echo ; echo
      ((jobindex++))

    done< <(grep "/000${runNumber}/" ${inputFileList})
    
    #CPass0 list of Calib files to merge
    #arr_cpass0_calib_list[${runNumber}]="${commonOutputPath}/meta/cpass0.calib.run${runNumber}.list"
    arr_cpass0_calib_list[${runNumber}]="cpass0.calib.run${runNumber}.list"
    echo "### Produces the list of CPass0 files to merge (executes locally) ###"
    echo "${arr_cpass0_calib_list[${runNumber}]}: benchmark.sh ${arr_cpass0_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues calibfile ${arr_cpass0_calib_list[${runNumber}]} ${arr_cpass0_outputs[*]}"
    echo ; echo

    #CPass0 merging
    echo "### Merges CPass0 files ###"
    #arr_cpass0_merged[${runNumber}]="${commonOutputPath}/meta/merge.cpass0.run${runNumber}.done"
    arr_cpass0_merged[${runNumber}]="merge.cpass0.run${runNumber}.done"
    echo "${arr_cpass0_merged[${runNumber}]}: benchmark.sh ${configFile} ${arr_cpass0_calib_list[${runNumber}]} ${copyFiles[@]}"
    echo -e "\t${alirootEnv} ./benchmark.sh MergeCPass0 \$OUTPATH/000${runNumber}/cpass0 ${currentDefaultOCDB} ${configFile} ${runNumber} ${arr_cpass0_calib_list[${runNumber}]} ${extraOpts[@]}"" "
    echo ; echo

    #CPass1 list of Calib/QA/ESD/filtered files
    # the trick with QA is to have the string "Stage.txt" in the file name of the list of directories with QA output to trigger
    # the production of the QA trending tree (only then the task->Finish() will be called in QAtrain_duo.C, on the grid
    # this corresponds to the last merging stage)
    #arr_cpass1_QA_list[${runNumber}]="${commonOutputPath}/meta/cpass1.QA.run${runNumber}.lastMergingStage.txt.list"
    arr_cpass1_QA_list[${runNumber}]="cpass1.QA.run${runNumber}.lastMergingStage.txt.list"
    echo "### Lists CPass1 QA ###"
    echo "${arr_cpass1_QA_list[${runNumber}]}: benchmark.sh ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues dir ${arr_cpass1_QA_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    #arr_cpass1_calib_list[${runNumber}]="${commonOutputPath}/meta/cpass1.calib.run${runNumber}.list"
    arr_cpass1_calib_list[${runNumber}]="cpass1.calib.run${runNumber}.list"
    echo "### Lists CPass1 Calib ###"
    echo "${arr_cpass1_calib_list[${runNumber}]}: benchmark.sh ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues calibfile ${arr_cpass1_calib_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo


  
    #CPass1 merging
    #arr_cpass1_merged[${runNumber}]="${commonOutputPath}/meta/merge.cpass1.run${runNumber}.done"
    arr_cpass1_merged[${runNumber}]="merge.cpass1.run${runNumber}.done"
    echo "### Merges CPass1 files ###"
    echo "${arr_cpass1_merged[${runNumber}]}: benchmark.sh ${configFile} ${arr_cpass1_calib_list[${runNumber}]} ${arr_cpass1_QA_list[${runNumber}]} ${copyFiles[@]}"
    echo -e "\t${alirootEnv} ./benchmark.sh MergeCPass1 \$OUTPATH/000${runNumber}/cpass1 ${currentDefaultOCDB} ${configFile} ${runNumber} ${arr_cpass1_calib_list[${runNumber}]} ${arr_cpass1_QA_list[${runNumber}]} ${arr_cpass1_filtered_list[${runNumber}]} ${extraOpts[@]}"" "
    echo ; echo

    #CPass0 wrapped in a profiling tool (valgrind,....)
    if [[ -n ${profilingCommand} ]]; then
      inputFile=$(grep -m1 "${runNumber}/" ${inputFileList})
      [[ -z ${nEventsProfiling} ]] && nEventsProfiling=2
      currentDefaultOCDB=$(setYear ${inputFile} ${defaultOCDB})
      jobindex="profiling"

      #arr_cpass0_profiled_outputs[${runNumber}]="${commonOutputPath}/meta/profiling.cpass0.job${jobindex}.run${runNumber}.done"
      arr_cpass0_profiled_outputs[${runNumber}]="profiling.cpass0.job${jobindex}.run${runNumber}.done"
      echo "### CPass0 in a profiler ###"
      echo "${arr_cpass0_profiled_outputs[${runNumber}]}: benchmark.sh ${configFile} ${copyFiles[@]}"
      profilingCommand=$(encSpaces "${profilingCommand}")
      echo -e "\t${alirootEnv} ./benchmark.sh CPass0 \$OUTPATH/000${runNumber}/${jobindex} ${inputFile} ${nEventsProfiling} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]} useProfilingCommand=${profilingCommand}"
      echo ; echo
    fi

  done #runs

  #Summary
  echo "### Summary ###"
  echo "summary.log: benchmark.sh ${configFile} ${arr_cpass0_outputs[*]} ${arr_cpass0_merged[*]} ${arr_cpass1_outputs[*]} ${arr_cpass1_merged[*]}"
  echo -e "\t${alirootEnv} ./benchmark.sh MakeSummary ${configFile} ${extraOpts[@]}"
  echo ; echo

  return 0
)

goPrintValues()
(
  #print the values given the key from any number of files (store in output file on second argument)
  if [[ $# -lt 3 ]]; then
    echo "goPrintValues key outputFile inputFiles"
    echo "if outputFile is \"-\" don't redirect to a file"
    return
  fi
  key=${1}
  outputFile=${2}
  [[ ${outputFile} == "-" ]] && outputFile=""
  shift 2
  cat $(parseListOfFiles "$@") | awk -v key=${key} '$0 ~ key" " {print $2}' | tee ${outputFile}
  return 0
)

goCreateQAplots()
(
  umask 0002
  mergedQAfileList=${1}
  productionID=${2}
  outputDir=${3}
  configFile=${4}
  shift 4
  alilog_info  "[BEGIN] goCreateQAplots() with following parameters $*"
  echo "$@"
  if ! parseConfig configFile=${configFile} "$@"; then 
     alilog_error "goCreateQAplots() Parsing config error [Paremeters] $*"
     return 1; 
  fi
  
  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  [[ -f ${alirootSource} && -z ${ALICE_ROOT} ]] && source ${alirootSource}

  [[ -z ${logOutputDir} ]] && logOutputDir=${PWD}
  [[ -z ${dontRedirectStdOutToLog} ]] && exec &> ${logOutputDir}/makeQAplots.log
  echo "${0} $*"

  mkdir -p ${outputDir}
  cd ${outputDir}
  [[ ! "${PWD}" =~ "${outputDir}" ]] && echo "PWD is not equal to outputDir=${outputDir}" && cd ${batchWorkingDirectory} && return 1

  inputFiles=(
              "${batchWorkingDirectory}/runQA.sh"
              "${ALICE_PHYSICS}/PWGPP/QA/scripts/runQA.sh"
  )
  for file in ${inputFiles[*]}; do
    [[ ! -f ${file##*/} && -f ${file} ]] && echo "copying ${file}" && cp -f ${file} .
  done

  echo "running QA with command:"
  echo ./runQA.sh inputList="${mergedQAfileList}" inputListHighPtTrees="${filteringList}" ocdbStorage="${defaultOCDB}"
  ./runQA.sh inputList="${mergedQAfileList}" inputListHighPtTrees="${filteringList}" ocdbStorage="${defaultOCDB}"
  cd ${batchWorkingDirectory}
  alilog_info  "[END] goCreateQAplots() with folloing parameters: $*"
  echo "$@"
  return 0
)

goTest()
(
  echo AA
)

alirootInfo()
(
  umask 0002
  # save aliroot repository info
  [[ -z "${ALICE_ROOT}" ]] && return 1
  
  echo "\${ALICE_ROOT}=${ALICE_ROOT}"
  echo "\${ROOTSYS}=${ROOTSYS}"
  echo "\${PATH}=${PATH}"
  echo "\${LD_LIBRARY_PATH}=${LD_LIBRARY_PATH}"
  echo
  
  pushd ${PWD}
  cd ${ALICE_ROOT}

  currentBranch=$(git rev-parse --abbrev-ref HEAD)
  git status
  echo ""
  echo ""
  git diff ${currentBranch}
  popd
  return 0
)

spitOutLocalOCDBaccessConfig()
{
  umask 0002
  #find ${1} -name "*root" | \
  /bin/ls -1 ${1}/*/*/*/*.root 2>/dev/null | \
  while read line
  do 
    local tmp=${line#${1}}
    echo ${tmp%/*} | \
    awk -v ocdb=${1} '{print "  man->SetSpecificStorage(\""$1"\",\"local://"ocdb"\");"}'
  done
  return 0
}

goMakeLocalOCDBaccessConfig()
{
  umask 0002
  # make a script that sets the specific storages form all the root files produced by CPass0
  local localOCDBpathCPass0=${1}
  local OCDBpathPrefix=${2-.}

  if [[ -f ${localOCDBpathCPass0} && ${localOCDBpathCPass0} =~ \.tgz$ ]]; then
    tar xzvvf ${localOCDBpathCPass0}
    local localOCDBpathCPass0="${OCDBpathPrefix}/OCDB"
  fi

  echo
  echo creating the specific storage script
  echo   localOCDBaccessConfig.C
  echo   based on OCDB: ${localOCDBaccessConfig}
  echo

  local tempLocalOCDB=""
  if [[ -f localOCDBaccessConfig.C ]]; then
    tempLocalOCDB=$(mktemp -t tempLocalOCDB.XXXXXX)
    egrep "SetSpecificStorage" localOCDBaccessConfig.C > ${tempLocalOCDB}
  fi

  echo "localOCDBaccessConfig()"                               >  localOCDBaccessConfig.C
  echo "{"                                                     >> localOCDBaccessConfig.C
  echo "  AliCDBManager* man = AliCDBManager::Instance();"     >> localOCDBaccessConfig.C
  spitOutLocalOCDBaccessConfig ${localOCDBpathCPass0}|sort|uniq  >> localOCDBaccessConfig.C
  [[ -f "${tempLocalOCDB}" ]] && cat ${tempLocalOCDB}              >> localOCDBaccessConfig.C
  echo "}"                                                     >> localOCDBaccessConfig.C

  [[ -f "${tempLocalOCDB}" ]] && rm -f ${tempLocalOCDB}

  if ! grep SetSpecificStorage localOCDBaccessConfig.C; then 
    echo
    echo "!!!!!!! CPass0 produced no OCDB entries"
    return 1
  fi
  return 0
}

goMakeFilteredTrees()
(
  outputDir=${1}
  runNumber=${2}
  #get path to input list
  inputListfiles=${3}
  #get scale number for tracks
  filterT=${4}
  #get scale number for V0s
  filterV=${5}
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
  alilog_info  "[BEGIN] goMakeFilteredTrees() with following parameters $*"
  shift 12
  if ! parseConfig configFile=${configFile} "$@"; then return 1; fi
  
  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  [[ -z ${commonOutputPath} ]] && commonOutputPath=${PWD}
  doneFileBase=filtering.cpass1.run${runNumber}.done
  doneFileTmp=${batchWorkingDirectory}/${doneFileBase}
  doneFile=${commonOutputPath}/meta/${doneFileBase}

  cat > filtering.log << EOF
  goMakeFilteredTrees config:
  runpath=${runpath}
  outputDir=${outputDir}
  commonOutputPath=${commonOutputPath}
  ALICE_ROOT=${ALICE_ROOT}
  PATH=${PATH}
  offsetEvent=$offsetEvent
  configFile=$configFile
  esdFileName=$esdFileName
  inputListfiles=$inputListfiles
  doneFile=$doneFile
EOF

  #runpath=${outputDir}
  #[[ ${reconstructInTemporaryDir} -eq 1 && -n ${TMPDIR} ]] && runpath=${TMPDIR}
  #[[ ${reconstructInTemporaryDir} -eq 1 ]] && runpath=$(mktemp -d -t goMakeFilteredTrees.XXXXXX)
  #mkdir -p ${outputDir}
  #mkdir -p ${runpath}
  #if ! cd ${runpath}; then 
  #  echo "PWD=$PWD is not the runpath=${runpath}"
  #  touch ${doneFile}
  #  return 1
  #fi
  
  if [[ -z ${pretend} ]];then
    aliroot -l -b -q "${ALICE_PHYSICS}/PWGPP/macros/runFilteringTask.C(\"${inputListfiles}\",${filterT},${filterV},\"${OCDBpath}\",${maxFiles},${offsetFile},${maxEvents},${offsetEvent},\"${esdFileName}\")" &>> filtering.log
  else
    sleep ${pretendDelay}
    touch filtering.log FilterEvents_Trees.root
  fi
  pwd
  /bin/ls
  summarizeLogs >>  ${doneFile}
  
  #echo mv -f * ${outputDir}
  #mv -f * ${outputDir}
  #[[ -f ${outputDir}/FilterEvents_Trees.root ]] && echo "filteredTree ${outputDir}/FilterEvents_Trees.root" >> ${doneFile}
  #cd ${commonOutputPath}
  #[[ "${runpath}" != "${outputDir}" ]] && rm -rf ${runpath}
 
  cp "$doneFileTmp" "$doneFile" || rm -f "$doneFileTmp" "$doneFile"
  [[ -n ${removeTMPdoneFile} ]] && rm -f ${doneFileTmp}
  alilog_info  "[END] goMakeFilteredTrees() with following parameters $*"
  return 0
)

submit()
{
  umask 0002
  [[ $# -lt 5 ]] && echo "at least 5 args needed, you supplied $#" && return 1
  JobID=${1}
  startID=${2}
  endID=${3}
  waitForJOBID=${4}
  command=${5}
  shift 5
  local commandArgs=("$@")

  #add quote strings around the extra arguments
  for ((i=0; i<${#commandArgs[@]}; i++)); do 
    commandArgs[i]=\"${commandArgs[i]}\"
  done

  [[ -z ${waitForJOBID} ]] && waitForJOBID=0

  newFarm=$(which qsub|grep "^/usr/bin/qsub")
  
  batchSystem="SGE"

  if [[ -z "${newFarm}" ]]
  then
    #old LSF
    # submit it (as job array)
    nFiles=$(( ${endID}-${startID}+1 ))
    while [ ${startID} -le ${nFiles}  ] ; do
      if [ $(expr ${nFiles} - ${startID}) -gt 999 ] ; then 
        endID=$(expr ${startID} + 999)
      else
        endID=${nFiles}
      fi      
    if [[ ${waitForJOBID} -eq 0 ]]; then
        echo ${batchCommand} ${batchFlags} -J "${JobID}[${startID}-${endID}]" -e "${commonOutputPath}/logs/job_%I.err" -o "${commonOutputPath}/logs/job_%I.out" "${command}"     
        ${batchCommand} ${batchFlags} -J "${JobID}[${startID}-${endID}]" -e "${commonOutputPath}/logs/job_%I.err" -o "${commonOutputPath}/logs/job_%I.out" "${command}"     
      else
        echo ${batchCommand} ${batchFlags} -J "${JobID}[${startID}-${endID}]" -w "ended(${waitForJOBID})" -e "${commonOutputPath}/logs/job_%I.err" -o "${commonOutputPath}/logs/job_%I.out" "${command}"     
        ${batchCommand} ${batchFlags} -J "${JobID}[${startID}-${endID}]" -w "ended(${waitForJOBID})" -e "${commonOutputPath}/logs/job_%I.err" -o "${commonOutputPath}/logs/job_%I.out" "${command}"     
      fi
      startID=$(expr ${endID} + 1)
    done
  else 
    #new SGE farm
    if [[ ${waitForJOBID} =~ "000" ]]; then
      echo ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
      ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
    else
      echo ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -hold_jid "${waitForJOBID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
      ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -hold_jid "${waitForJOBID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
    fi
  fi
  return 0
}

goSubmitBatch()
{
  if [[ $# -lt 3 ]]; then
    echo "minimal use:"
    echo " ${0} submit fileList productionID configFile"
    return 0
  fi

  productionID=${1}
  inputList=${2}
  configFile=${3:-"benchmark.config"}
  #if which greadlink; then configFile=$(greadlink -f ${configFile}); fi
  shift 3
  extraOpts=("$@")
  if ! parseConfig configFile=${configFile} "${extraOpts[@]}"; then return 1; fi
  
  #batch systems/makeflow sometimes handle spaces in arguments poorly, so encode them
  for (( i=0;i<${#extraOpts[@]};i++ )); do 
    extraOpts[i]=$(encSpaces "${extraOpts[i]}")
  done
  extraOpts+=("encodedSpaces=1")
  #this removes the copy of the done file used by makeflow (in the running dir)
  extraOpts+=("removeTMPdoneFile=1")

  #record the working directory provided by the batch system
  batchWorkingDirectory=${PWD}

  #redirect all output to submit.log
  echo "redirecting all output to ${PWD}/submit_${productionID//"/"/_}.log"
  exec 7>&1
  exec 1>submit_${productionID//"/"/_}.log 2>&1

  umask 0002
  echo ${0}" submit $*"
  if [[ -z "${inputList}" || -z "${productionID}" ]]
  then
    echo
    echo " Usage: ${0} submit inputList productionID [configFile=benchmark.config]"
    echo
    return
  fi

  # check if config file is there
  if [ ! -f ${configFile} ]; then
    echo "ERROR! Config File '${configFile}' not found" >&2
    return
  else
    echo "Using Config File: '${configFile}'"
  fi

  [[ ! -f ${alirootEnv} ]] && echo "alirootEnv script ${alirootEnv} not found!..." && return 1

  #move the script, config and some other stuff to ${commonOutputPath} first, then use them from there
  self=${0}
  #if which greadlink; then self=$(greadlink -f "${0}"); fi
  configPath=$(dirname ${configFile})
  export commonOutputPath=${baseOutputDirectory}/${productionID}
  
  mkdir -p ${commonOutputPath}
  mkdir -p ${commonOutputPath}/logs
  mkdir -p ${commonOutputPath}/meta

  cp ${self} ${commonOutputPath}
  cp ${configFile} ${commonOutputPath}
  cp ${inputList} ${commonOutputPath}
  self=${commonOutputPath}/${self##*/}
  chmod u+x ${self}
  configFile=${commonOutputPath}/${configFile##*/}
  inputList=${commonOutputPath}/${inputList##*/}

  #convert to absolut pathnames
  #if which greadlink; then inputList=$(greadlink -f "${inputList}"); fi
  #make list of runs
  if [[ -z ${runNumber} ]]; then
    listOfRuns=($(while read x; do guessRunNumber ${x}; done < ${inputList} | sort | uniq))
  else
    listOfRuns=${runNumber}
  fi

  #if which greadlink; then alirootSource=$(greadlink -f "${alirootSource}"); fi

  echo ""
  echo "### BEGIN CONFIGURATION ###"
  echo ""
  echo "GENERAL:"
  echo ""
  echo "    productionID:    ${productionID}"
  echo "    batchCommand:    ${batchCommand}"
  echo "    batchFlags:      ${batchFlags}"
  echo "    alirootEnv:   ${alirootEnv}"
  ${alirootEnv} echo '    ALICE_ROOT:      ${ALICE_ROOT}'
  ${alirootEnv} echo '    ALIROOT_RELEASE: ${ALICE_RELEASE}'
  echo "    inputList:       ${inputList}"
  echo "    configPath:      ${configPath}"
  echo "    commonOutputPath:      ${commonOutputPath}"
  echo "    defaultOCDB:     ${defaultOCDB}"
  echo "      autoOCDB: ${autoOCDB}"
  echo "    recoTriggerOptions:   ${recoTriggerOptions}"
  echo "    runs:"
  echo "      ${listOfRuns[*]}"
  echo ""
  echo "THE TRAIN WILL RUN:"

  if [ ${runCPass0reco} -eq 1 ]; then
    echo "    Pass0 - Recontruction"
  fi

  if [ ${runCPass0MergeMakeOCDB} -eq 1 ]; then
    echo "    Pass0 - merging and OCDB export"
  fi

  if [ ${runCPass1reco} -eq 1 ]; then
    echo "    Pass1 - Recontruction"
  fi
  if [ ${runCPass1MergeMakeOCDB} -eq 1 ]; then
    echo "    Pass1 - merging and OCDB export"
  fi

  echo ""
  echo "LIMITS:"
  echo "    max. Events/Chunk:   ${nEvents}"
  echo "    max. Number of Chunks per Run:     ${nMaxChunks}"
  echo ""
  echo "### END CONFIGURATION ###"
  echo ""


  # check if input file is there
  if [ ! -f ${inputList} ]; then
    echo "ERROR! Input List '${inputList}' not found" >&2
    return
  fi

  # define jobid (for dependent jobs)
  date=$(date +%Y_%m_%d_%H%M%S)
  #for each run we submit one jobarray:
  for runNumber in ${listOfRuns[*]}; do
    
    [[ -z ${runNumber} ]] && continue
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

    oneInputFile=$(egrep -m1 "${runNumber}/" ${inputList})

    currentDefaultOCDB=${defaultOCDB}
    [[ -z ${autoOCDB} ]] && autoOCDB=1
    if [[ ${autoOCDB} -ne 0 ]]; then
      currentDefaultOCDB=$(setYear ${oneInputFile} ${defaultOCDB})
    fi
    period=$(guessPeriod ${oneInputFile})
    year=$(guessYear ${oneInputFile})

    echo "submitting run ${runNumber} with OCDB ${currentDefaultOCDB}"

    ###############################################################################
    #run one chunk with valgrind:
    if [[ -n ${profilingCommand} ]]; then
      [[ -z ${nEventsProfiling} ]] && nEventsProfiling=2
      [[ -z ${profilingCommand} ]] && profilingCommand="/usr/bin/valgrind --tool=callgrind --num-callers=40 -v --trace-children=yes"
      submit "profile-${JOBpostfix}" 1 1 000 "${alirootEnv} ${self}" CPass0 ${commonOutputPath}/${year}/${period}/000${runNumber}/${jobindex} ${oneInputFile} ${nEventsProfiling} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} useProfilingCommand=$(encSpaces "${profilingCommand}") "${extraOpts[@]}"
    fi 

    ################################################################################
    ################################################################################
    # run the CPass0 if requested

    if [ ${runCPass0reco} -eq 1 ]; then

      echo
      echo "starting CPass0... for run ${runNumber}"
      echo

      # create directory and copy all files that are needed
      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass0"
      mkdir -p ${targetDirectory}

      filesCPass0=( 
                    "${configPath}/runCPass0.sh"
                    "${configPath}/recCPass0.C"
                    "${configPath}/runCalibTrain.C"
                    "${configPath}/localOCDBaccessConfig.C"
                    "${configPath}/OCDB*.root"
                    "${configPath}/sim.C"
                    "${configPath}/rec.C"
                    "${configPath}/Config.C"
      )
      for file in ${filesCPass0[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      localInputList=${targetDirectory}/${inputList##*/}
      [[ ! -f ${localInputList} ]] && egrep "\/000${runNumber}\/" ${inputList} > ${localInputList}
      # limit nFiles to nMaxChunks
      nFiles=$(wc -l < ${localInputList})
      [[ ${nFiles} -eq 0 ]] && echo "list contains ZERO files! exiting..." && return 1
      echo "raw files in list:    ${nFiles}"
      if [[ ${nMaxChunks} -gt 0 && ${nMaxChunks} -le ${nFiles} ]]; then
        nFiles=${nMaxChunks}
      fi
      echo "raw files to process: ${nFiles}"
      [[ -z "${percentProcessedFilesToContinue}" ]] && percentProcessedFilesToContinue=100
      if [[ ${percentProcessedFilesToContinue} -eq 100 ]]; then
        nFilesToWaitFor=${nFiles}
      else
        nFilesToWaitFor=$(( ${nFiles}-${nFiles}/(100/(100-${percentProcessedFilesToContinue})) ))
      fi
      echo "requested success rate is ${percentProcessedFilesToContinue}%"
      echo "merging will start after ${nFilesToWaitFor} jobs are done"

      submit ${JOBID1} 1 ${nFiles} 000 "${alirootEnv} ${self}" CPass0 ${targetDirectory} ${localInputList} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} -1 "${extraOpts[@]}"

      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "${JOBID1wait}" 1 1 000 "${alirootEnv} ${self}" WaitForOutput ${commonOutputPath} "meta/cpass0.job\*.run${runNumber}.done" ${nFilesToWaitFor} ${maxSecondsToWait}
      LASTJOB=${JOBID1wait}

    fi #end running CPass0
    ################################################################################


    ################################################################################
    # submit merging of CPass0, depends on the reconstruction

    if [ ${runCPass0MergeMakeOCDB} -eq 1 ]; then

      echo
      echo "submit CPass0 merging for run ${runNumber}"
      echo

      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass0"
      mkdir -p ${targetDirectory}

      #copy the scripts
      filesMergeCPass0=(
                        "${configPath}/OCDB.root"
                        "${configPath}/mergeMakeOCDB.byComponent.sh"
                        "${configPath}/mergeMakeOCDB.sh"
                        "${configPath}/localOCDBaccessConfig.C"
                        "${configPath}/mergeByComponent.C"
                        "${configPath}/makeOCDB.C"
                        "${configPath}/merge.C"
      )
      for file in ${filesMergeCPass0[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done
  
      submit ${JOBID2} 1 1 "${LASTJOB}" "${alirootEnv} ${self}" MergeCPass0 ${targetDirectory} ${currentDefaultOCDB} ${configFile} ${runNumber} cpass0.calib.run${runNumber}.list "${extraOpts[@]}"
      LASTJOB=${JOBID2}

      if [[ -n ${generateMC} ]]; then
        submit "mrl${JOBpostfix}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" PrintValues sim ${commonOutputPath}/meta/sim.run${runNumber}.list ${commonOutputPath}/meta/cpass0.job*.run${runNumber}.done
        LASTJOB="mrl${JOBpostfix}"
      fi

      echo
    fi
    # end of merging CPass0
    ################################################################################

    ################################################################################
    ################################################################################
    # run the CPass1 if requested

    if [ ${runCPass1reco} -eq 1 ]; then

      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass1"
      rm -f ${commonOutputPath}/meta/cpass1.job*.run${runNumber}.done

      # safety feature: if we are re-running for any reason we want to delete the previous output first.
      [[ -d ${targetDirectory} ]] && rm -rf ${targetDirectory}/* && echo "removed old output at ${targetDirectory}/*"

      echo
      echo "starting CPass1... for run ${runNumber}"
      echo

      # create directory and copy all files that are needed
      mkdir -p ${targetDirectory}
      
      filesCPass1=( 
                    "${configPath}/runCPass1.sh"
                    "${configPath}/recCPass1.C"
                    "${configPath}/recCPass1_OuterDet.C"
                    "${configPath}/runCalibTrain.C"
                    "${configPath}/QAtrain_duo.C"
                    "${configPath}/mergeQAgroups.C"
                    "${configPath}/localOCDBaccessConfig.C"
                    "${configPath}/OCDB.root"
      )
      for file in ${filesCPass1[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      if [[ -n ${generateMC} ]]; then
        localInputList=${commonOutputPath}/meta/sim.run${runNumber}.list
      else
        localInputList=${targetDirectory}/${inputList##*/}
        [[ ! -f ${localInputList} ]] && egrep "\/000${runNumber}\/" ${inputList} > ${localInputList}
      fi
      # limit nFiles to nMaxChunks
      nFiles=$(wc -l < ${localInputList})
      [[ ${nFiles} -eq 0 ]] && echo "list contains ZERO files! continuing..." && continue
      echo "raw files in list:    ${nFiles}"
      if [[ ${nMaxChunks} -gt 0 && ${nMaxChunks} -le ${nFiles} ]]; then
        nFiles=${nMaxChunks}
      fi
      echo "raw files to process: ${nFiles}"
      [[ -z "${percentProcessedFilesToContinue}" ]] && percentProcessedFilesToContinue=100
      if [[ ${percentProcessedFilesToContinue} -eq 100 ]]; then
        nFilesToWaitFor=${nFiles}
      else
        nFilesToWaitFor=$(( ${nFiles}-${nFiles}/(100/(100-${percentProcessedFilesToContinue})) ))
      fi
      echo "requested success rate is ${percentProcessedFilesToContinue}%"
      echo "merging will start after ${nFilesToWaitFor} jobs are done"

      submit ${JOBID4} 1 ${nFiles} "${LASTJOB}" "${alirootEnv} ${self}" CPass1 ${targetDirectory} ${localInputList} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} -1 "${extraOpts[@]}"

      ################################################################################
      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "${JOBID4wait}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" WaitForOutput ${commonOutputPath} "meta/cpass1.job\*.run${runNumber}.done" ${nFilesToWaitFor} ${maxSecondsToWait}
      LASTJOB=${JOBID4wait}
      ################################################################################

      echo
    fi #end running CPass1

    ################################################################################
    # submit merging of CPass1, depends on the reconstruction
    if [ ${runCPass1MergeMakeOCDB} -eq 1 ]; then

      echo
      echo "submit CPass1 merging for run ${runNumber}"
      echo

      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass1"
      rm -f ${commonOutputPath}/meta/merge.cpass1.run${runNumber}.done
      mkdir -p ${targetDirectory}

      # copy files 
      filesMergeCPass1=(
                        "${configPath}/OCDB.root"
                        "${configPath}/localOCDBaccessConfig.C"
                        "${configPath}/mergeMakeOCDB.byComponent.sh"
                        "${configPath}/mergeByComponent.C"
                        "${configPath}/makeOCDB.C"
                        "${configPath}/merge.C"
                        "${configPath}/mergeMakeOCDB.sh"
                        "${configPath}/QAtrain_duo.C"
			"${configPath}/mergeQAgroups.C"
      )
      for file in ${filesMergeCPass1[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      submit "${JOBID5}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" MergeCPass1 ${targetDirectory} ${currentDefaultOCDB} ${configFile} ${runNumber} cpass1.calib.run${runNumber}.list cpass1.QA.run${runNumber}.lastMergingStage.txt.list cpass1.filtered.run${runNumber}.list "${extraOpts[@]}"
      LASTJOB=${JOBID5}
      echo
    fi

    ###############################
    #if [ ${runESDfiltering} -eq 1 ]; then
    #  rm -f ${commonOutputPath}/cpass1.ESD.run${runNumber}.list
    #  rm -f ${commonOutputPath}/meta/filtering.cpass1.run*.done
    #  echo
    #  echo submitting filtering for run ${runNumber}
    #  echo
    #  submit "${JOBmakeESDlistCPass1}" 1 1 "${LASTJOB}" "${self}" PrintValues esd ${commonOutputPath}/meta/cpass1.ESD.run${runNumber}.list ${commonOutputPath}/meta/cpass1.job*.run${runNumber}.done 
    #  submit "${JOBfilterESDcpass1}" 1 1 "${JOBmakeESDlistCPass1}" "${alirootEnv} ${self}" MakeFilteredTrees ${commonOutputPath}/${year}/${period}/000${runNumber}/cpass1 ${runNumber} ${commonOutputPath}/meta/cpass1.ESD.run${runNumber}.list ${filteringFactorHighPt} ${filteringFactorV0s} ${currentDefaultOCDB} 1000000 0 10000000 0 ${configFile} AliESDs_Barrel.root "${extraOpts[@]}"
    #  LASTJOB=${JOBfilterESDcpass1}
    #fi

  done

  #################################################################################
  #################################################################################
  #if [ ${runESDfiltering} -eq 1 ]; then
  #  submit "${JOBID5wait}" 1 1 "${LASTJOB}" "${self}" WaitForOutput ${commonOutputPath} "meta/filtering.cpass1.run/\*.done" "${#listOfRuns[@]}" ${maxSecondsToWait}
  #else
    submit "${JOBID5wait}" 1 1 "${LASTJOB}" "${self}" WaitForOutput ${commonOutputPath} "meta/merge.cpass1.run\*.done" ${#listOfRuns[@]} ${maxSecondsToWait}
  #fi
  LASTJOB=${JOBID5wait}

  #################################################################################
  echo
  echo "submit make a summary"
  echo

  [[ -z ${alirootEnvQA} ]] && alirootEnvQA=$(encSpaces "${alirootEnv}")
  submit "${JOBID6}" 1 1 "${LASTJOB}" "${alirootEnvQA} ${self}" MakeSummary ${configFile} "commonOutputPath=${commonOutputPath}"
  LASTJOB=${JOBID6}
  #################################################################################
  
  #restore stdout
  exec 1>&7 7>&-
  echo "jobs submitted."
  return 0
}

goWaitForOutput()
(
  #
  # will be nice to make documentation
  # Used to fomaly define dependencies of jobs (wait until something finish) 
  # Based on my experience (MI) with benchmark - we are waiting for particualr files to exist or timeout
  # In case of infinite waiting time (because of some misbehaving) easiest way to process is to kill job executing this command
  #      (jobs on batch farm has PID dependency)
  #
  
  # Action:
  # Input:
  # Output:

  umask 0002
  [[ $# -lt 3 ]] && echo "goWaitForOutput() wrong number of arguments, exiting.." && return 1
  echo Start:goWaitForOutput
  echo searchPath=${1}
  echo fileName=${2}
  echo numberOfFiles=${3}
  echo maxSecondsToWait=${4}
  searchPath=${1}
  fileName=${2}
  numberOfFiles=${3}
  maxSecondsToWait=${4}
  echo "command to be executed: /bin/ls -1 ${searchPath}/${fileName}"
  [[ -z "${maxSecondsToWait}" ]] && maxSecondsToWait=$(( 3600*12 ))
  while true; do
    n=$(/bin/ls -1 ${searchPath/\\/}/${fileName/\\/} 2>/dev/null | wc -l)
    [[ ${n} -gt 0 ]] && echo "found ${n} X ${fileName}"
    [[ ${n} -ge ${numberOfFiles} ]] && break
    [[ ${SECONDS} -gt ${maxSecondsToWait} ]] && echo "timeout of ${maxSecondsToWait}s!" && break
    sleep 60
  done
  echo "DONE! exiting..."
  echo End:goWaitForOutput
  return 0
)

goMakeMergedSummaryTree()
(
  # create list of calibration entries
  # takes no arguments, just run it in the base output
  # directory with the following files in the working directory
  # 
  #  Calibration file lists:
  #       cpass0.dcsTree.list, cpass1.dcsTree.list
  #  QA trending root file:
  #       trending.root
  # 
  #  Production infoo ascii files:
  #       summary_pass0.tree
  #       summary_pass1.tree
  #    
    alilog_info  "[BEGIN] goMakeMergedSummaryTree() with following parameters $*"
  [[ ! -f cpass0.dcsTree.list ]] && echo "no cpass0.dcsTree.list" && return 1
  [[ ! -f cpass1.dcsTree.list ]] && echo "no cpass1.dcsTree.list" && return 1
  [[ ! -f trending.root ]] && echo "no trending.root" && return 1
  [[ ! -f summary_pass0.tree ]] && echo "no summary_pass0.tree" && return 1
  [[ ! -f summary_pass1.tree ]] && echo "no summary_pass1.tree" && return 1

  #first, dump the C macro to file
  cat << EOF > mergeTree.C
  //
  // Merge summary information
  // Following files are expected to be in the working directory
  //
  // Calibration file lists:
  //      cpass0.dcsTree.list, cpass1.dcsTree.list
  // QA trending root file:
  //      trending.root
  //
  // Production infoo ascii files:
  //      summary_pass0.tree
  //      summary_pass1.tree
  //   
  void mergeTree(){
    //
    //
    //
    // Calibration values dump
    //
    //Printf("MakeTreeFromList cpass0.dcsTree.list");
    AliXRDPROOFtoolkit::MakeTreeFromList("Calib.TPC.CPass0.root", "dcs","dcs","cpass0.dcsTree.list",1);
    //Printf("MakeTreeFromList cpass1.dcsTree.list");
    AliXRDPROOFtoolkit::MakeTreeFromList("Calib.TPC.CPass1.root", "dcs","dcs","cpass1.dcsTree.list",1);
    //
    // Calibration status dump
    //
    TFile *fprod = TFile::Open("fproduction.root","recreate");
    TTree  tree0, tree1;
    //Printf("reading summary_pass0.tree");
    tree0.ReadFile("summary_pass0.tree");
    //Printf("reading summary_pass1.tree");
    tree1.ReadFile("summary_pass1.tree");
    tree0.Write("CPass0");
    tree1.Write("CPass1");
    fprod->Close();
    //
    //
    //
    TString stringSetup="";
    stringSetup+="1#QA.TPC#run#SummaryTPCQA/tpcQA#trending.root+";  // 
    stringSetup+="1#Calib.TPC.CPass0#run#dcs#Calib.TPC.CPass0.root+";  // 
    stringSetup+="1#Calib.TPC.CPass1#run#dcs#Calib.TPC.CPass1.root+";  // 
    //
    stringSetup+="1#CPass0#runnumber#CPass0#fproduction.root+";  // 
    stringSetup+="1#CPass1#runnumber#CPass1#fproduction.root+";  // 
    //
    //Printf("stringSetup: %s", stringSetup.Data());
    AliXRDPROOFtoolkit::JoinTreesIndex("outAll.root","joinAll","run",stringSetup.Data(), 1);
  }
EOF

  aliroot -b -q "mergeTree.C" > mergeTrees.log
  alilog_info  "[END] goMakeMergedSummaryTree() with following parameters $*"
  return $?
)

goMakeSummary()
(
  # All the final stuff goes in here for ease of use:
  #  - summary logs
  #  - qa plot making
  #  - final file lists
  # Runs in current dir - in makeflow mode it can run LOCAL, then the QA plots and summaries
  # will appear in the submission dir.

  #some defaults:
  log=summary.log
  jsonLog=summary.json
  productionID=qa
  alilog_info "[BEGIN] goMakeSummary() with following parameters $*"
  configFile=$1
  shift 1
  extraOpts=("$@")

  parseConfig configFile=$configFile "${extraOpts[@]}" || return 1 
  batchWorkingDirectory=$PWD

  logTmp="${batchWorkingDirectory}/${log}"
  jsonLogTmp="${batchWorkingDirectory}/${jsonLog}"
  
  [[ -f "$alirootSource" && -z "$ALICE_ROOT" ]] && source "$alirootSource"
  [[ ! -f "$configFile" ]] && echo "no config file ${configFile}!" && return 1

  # If we use workqueue all the files are in PWD, on batch systems they are in the meta subdir.
  [[ -z ${commonOutputPath} ]] && commonOutputPath=${batchWorkingDirectory}
  metadir=${commonOutputPath}/meta
  [[ ! -d $metadir ]] && metadir="$PWD"

  exec &> >(tee ${logTmp})

  # Summarize the global stuff
  echo "env script: ${alirootSource} ${alirootEnv}"
  echo "ALICE_ROOT=${ALICE_ROOT}"
  echo "commonOutputPath=${commonOutputPath}"

  # json header: open array of objects
  echo '[' > "$jsonLogTmp"

  echo "total numbers for the production:"
  echo
  awk 'BEGIN {nFiles=0;nCore=0;} 
  /^calibfile/ {nFiles++;} 
  /core dumped/ {nCore++i;}
  END {print     "cpass0 produced "nFiles" calib files, "nCore" core files";}' "$metadir"/cpass0.job*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /\/rec.log OK/ {nOK++;} 
  /\/rec.log BAD/ {nBAD++;} 
  /stderr BAD/ {if ($0 ~ /rec.log/){nBAD++;}}
  END {print     "cpass0 reco:  OK: "nOK"\tBAD: "nBAD;}' "$metadir"/cpass0.job*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /\/calib.log OK/ {nOK++;} 
  /\/calib.log BAD/ {nBAD++;} 
  END {print "cpass0 calib: OK: "nOK"\tBAD: "nBAD;}' "$metadir"/cpass0.job*done 2>/dev/null

  awk 'BEGIN {nOK=0; nBAD=0; } 
  /merge.log OK/ {nOK++;} 
  /merge.log BAD/ {nBAD++;} 
  END {print "cpass0 merge: OK: "nOK"\tBAD: "nBAD;}' "$metadir"/merge.cpass0*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /ocdb.log OK/ {nOK++;} 
  /ocdb.log BAD/ {nBAD++;} 
  END {print   "cpass0 OCDB:  OK: "nOK"\tBAD: "nBAD;}' "$metadir"/merge.cpass0*done 2>/dev/null

  echo
  awk 'BEGIN {nFiles=0;nCore=0;} 
  /^calibfile/ {nFiles++;} 
  /core dumped/ {nCore++;}
  END {print     "cpass1 produced "nFiles" calib files, "nCore" core files";}' "$metadir"/cpass1.job*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /\/rec.log OK/ {nOK++;} 
  /\/rec.log BAD/ {nBAD++;} 
  /stderr BAD/ {if ($0 ~ /rec.log/){nBAD++;}}
  END {print     "cpass1 reco:  OK: "nOK"\tBAD: "nBAD;}' "$metadir"/cpass1.job*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /\/calib.log OK/ {nOK++;} 
  /\/calib.log BAD/ {nBAD++;} 
  END {print "cpass1 calib: OK: "nOK"\tBAD: "nBAD;}' "$metadir"/cpass1.job*done 2>/dev/null

  awk 'BEGIN {nOK=0; nBAD=0; } 
  /merge.log OK/ {nOK++;} 
  /merge.log BAD/ {nBAD++;} 
  END {print "cpass1 merge: OK: "nOK"\tBAD: "nBAD;}' "$metadir"/merge.cpass1*done 2>/dev/null
  awk 'BEGIN {nOK=0; nBAD=0; } 
  /ocdb.log OK/ {nOK++;} 
  /ocdb.log BAD/ {nBAD++;} 
  END {print   "cpass1 OCDB:  OK: "nOK"\tBAD: "nBAD;}' "$metadir"/merge.cpass1*done 2>/dev/null

  echo
  echo per run stats:
  /bin/ls -1 "$metadir"/merge.cpass0.run*.done | while read x; do

    dir=$(goPrintValues dir - ${x})
    runNumber=$(guessRunNumber ${dir})
    [[ -z ${runNumber} ]] && continue

    if $(/bin/ls "$metadir"/cpass0.job*.run${runNumber}.done &> /dev/null); then
      statusCPass0=( $(
      awk 'BEGIN {nOKrec=0;nBADrec=0;nOKcalib=0;nBADcalib=0;nOKstderr=0;nBADstderr=0;}
      /\/rec.log OK/ {nOKrec++;} 
      /\/rec.log BAD/ {nBADrec++;}
      /stderr BAD/ {if ($0 ~ /rec.log/) {nBADrec++;} nBADstderr++;}
      /stderr OK/ {nOKstderr++;}
      /\/calib.log OK/ {nOKcalib++;}
      /\/calib.log BAD/ {nBADcalib++}
      END {print ""nOKrec" "nBADrec" "nOKstderr" "nBADstderr" "nOKcalib" "nBADcalib;}' "$metadir"/cpass0.job*.run${runNumber}.done 2>/dev/null
      ) ) 
    fi

    if $(/bin/ls "$metadir"/cpass1.job*.run${runNumber}.done &>/dev/null); then
      statusCPass1=( $(
      awk 'BEGIN {nOKrec=0;nBADrec=0;nOKcalib=0;nBADcalib=0;nOKstderr=0;nBADstderr=0;nQAbarrelOK=0;nQAbarrelBAD=0;nQAouterOK=0;nQAouterBAD=0;}
      /\/rec.log OK/ {nOKrec++;} 
      /\/rec.log BAD/ {nBADrec++;}
      /stderr BAD/ {if ($0 ~ /rec.log/) nBADrec++;nBADstderr++;}
      /stderr OK/ {nOKstderr++;}
      /\/calib.log OK/ {nOKcalib++;}
      /\/calib.log BAD/ {nBADcalib++}
      /\/qa_barrel.log OK/ {nQAbarrelOK++;}
      /\/qa_barrel.log BAD/ {nQAbarrelBAD++;}
      /\/qa_outer.log OK/ {nQAouterOK++;}
      /\/qa_outer.log BAD/ {nQAouterBAD++;}
      END {print ""nOKrec" "nBADrec" "nOKstderr" "nBADstderr" "nOKcalib" "nBADcalib" "nQAbarrelOK" "nQAbarrelBAD" "nQAouterOK" "nQAouterBAD;}' "$metadir"/cpass1.job*.run${runNumber}.done 2>/dev/null
      ) ) 
    fi

    statusOCDBcpass0=$(awk '/ocdb.log/ {print $2} ' ${x} 2>/dev/null)
    statusOCDBcpass1=$(awk '/ocdb.log/ {print $2}' ${x/cpass0/cpass1} 2>/dev/null)
    statusQA=$(awk '/mergeMakeOCDB.log/ {print $2}' ${x/cpass0/cpass1} 2>/dev/null)

    printf "%s\t ocdb.log cpass0: %s\t ocdb.log cpass1: %s\tqa.log:%s\t| cpass0: rec:%s/%s stderr:%s/%s calib:%s/%s cpass1: rec:%s/%s stderr:%s/%s calib:%s/%s QAbarrel:%s/%s QAouter:%s/%s\n" \
           ${runNumber} ${statusOCDBcpass0} ${statusOCDBcpass1} ${statusQA} \
           ${statusCPass0[0]} ${statusCPass0[1]} ${statusCPass0[2]} ${statusCPass0[3]} ${statusCPass0[4]} ${statusCPass0[5]} \
           ${statusCPass1[0]} ${statusCPass1[1]} ${statusCPass1[2]} ${statusCPass1[3]} ${statusCPass1[4]} ${statusCPass1[5]} ${statusCPass1[6]} ${statusCPass1[7]} ${statusCPass1[8]} ${statusCPass1[9]}

    # produce json summary
    statusOCDBcpass0json=false
    statusOCDBcpass1json=false
    statusQAjson=false
    [[ "$statusOCDBcpass0" == 'OK' ]] && statusOCDBcpass0json=true
    [[ "$statusOCDBcpass1" == 'OK' ]] && statusOCDBcpass1json=true
    [[ "$statusQA" == 'OK' ]] && statusQAjson=true
    cat >> "$jsonLogTmp" <<EOF
    {
      run: ${runNumber},
      status: { ocdb_pass0: ${statusOCDBcpass0json}, ocdb_pass1: ${statusOCDBcpass1json}, qa: ${statusQAjson} },
      cpass0: {
        reco: { n_ok: ${statusCPass0[0]}, n_bad: ${statusCPass0[1]} },
        stderr: { n_ok: ${statusCPass0[2]}, n_bad: ${statusCPass0[3]} },
        calib: { n_ok: ${statusCPass0[4]}, n_bad: ${statusCPass0[5]} }
      },
      cpass1: {
        reco: { n_ok: ${statusCPass1[0]}, n_bad: ${statusCPass1[1]} },
        stderr: { n_ok: ${statusCPass1[2]}, n_bad: ${statusCPass1[3]} },
        calib: { n_ok: ${statusCPass1[4]}, n_bad: ${statusCPass1[5]} },
        qabarrel: { n_ok: ${statusCPass1[6]}, n_bad: ${statusCPass1[7]} },
        qarouter: { n_ok: ${statusCPass1[8]}, n_bad: ${statusCPass1[9]} }
      }
    },
EOF

  done

  # json footer: close array of objects
  echo ']' >> "${jsonLogTmp}"

  #make file lists - QA, trending, stacktraces, filtering and calibration
  ### wait for the merging of all runs to be over ###
  goPrintValues qafile remote.qa.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues calibfile remote.calib.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues trendingfile remote.trending.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues filteredTree remote.filtering.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues dcsTree remote.cpass0.dcsTree.list "$metadir"/merge.cpass0.run*.done &>/dev/null
  goPrintValues dcsTree remote.cpass1.dcsTree.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues stacktrace.log remote.cpass0.stacktrace.list "$metadir"/*cpass0*done &>/dev/null
  goPrintValues stacktrace.log remote.cpass1.stacktrace.list "$metadir"/*cpass1*done &>/dev/null
  goPrintValues esd remote.cpass1.esd.list "$metadir"/*cpass1*done &> /dev/null
 
  #copy all the files to a local dir tree
  for remoteList in remote.qa.list \
                    remote.calib.list \
                    remote.trending.list \
                    remote.filtering.list \
                    remote.cpass0.dcsTree.list \
                    remote.cpass1.dcsTree.list \
                    remote.cpass0.stacktrace.list \
                    remote.cpass1.stacktrace.list; do
    localList=${remoteList#remote.}
    rm -f "$localList" && touch "$localList"
    while read sourceFile; do
      destinationFile="${PWD}/${sourceFile#${baseOutputDirectory}}"
      copyFileFromRemote "$sourceFile" "$(dirname "${destinationFile}")" && \
        echo "$destinationFile" >> "$localList"
    done < <(cat "$remoteList")
  done

  #summarize the stacktraces
  stackTraceTree @stacktraces.cpass0.list > stacktrace_cpass0.tree
  stackTraceTree @stacktraces.cpass1.list > stacktrace_cpass1.tree

  #merge trending
  rm -f trending.root
  goMerge trending.list trending.root ${configFile} "${extraOpts[@]}" &> mergeTrending.log

  printExec goMakeSummaryTree "${metadir}" 0
  printExec goMakeSummaryTree "${metadir}" 1 

  # Checked until this point.

  goCreateQAplots "${PWD}/qa.list" "${productionID}" "QAplots" "${configFile}" "${extraOpts[@]}" filteringList="${PWD}/filtering.list" &>createQAplots.log

  #make a merged summary tree out of the QA trending, dcs trees and log summary trees
  goMakeMergedSummaryTree

  #if set, email the summary
  [[ -n ${MAILTO} ]] && cat ${logTmp} | mail -s "benchmark ${productionID} done" ${MAILTO}

  # Copy all, recursively.
  while read cpdir; do
    copyFileToRemote $cpdir/* $commonOutputPath/$cpdir
  done < <( find . -type d )

  alilog_info "[END] goMakeSummary() with following parameters $*"
  return 0
)

goMakeSummaryTree()
(
  alilog_info "[BEGIN] goMakeSummaryTree() with following parameters $*" 
  [[ $# -lt 1 ]] && return 1

  #1. define vars/arrays
  DIR=$1 #the directory with the meta files
  pass=${2-"0"} #pass from input
  outfile="summary_pass${pass}.tree"
  Ncolumns=0
  errfile=${outfile/tree/err}
  rm -f "$outfile" "$errfile"

  declare -a counterString=( TOFevents TOFtracks TPCevents TPCtracks TRDevents TRDtracks T0events \
                             SDDevents SDDtracks MeanVertexevents )
  Ncounter=${#counterString[@]}

  declare -a statusString=(TRDStatus TOFStatus TPCStatus T0Status MeanVertexStatus)
  Nstatus=${#statusString[@]}

  declare -a ratesString=(rec stderr calib qa_barrel qa_outer)
  Nrates=${#ratesString[@]}

  runs=( $(goPrintValues dir - merge.cpass0* | while read x; do guessRunNumber $x; done) )
  Nruns=${#runs[@]}

  echo -n runnumber/I >>${outfile}
  echo -n :cpass${pass}status/I >>${outfile}
  echo -n :cpass${pass}QAstatus/I >>${outfile}
  for i in ${ratesString[@]}; do
    echo -n :${i}OK/I >>${outfile}
    echo -n :${i}BAD/I >>${outfile}
  done
  for i in ${counterString[@]} ${statusString[@]} ; do
    echo -n :${i}/I >>${outfile}
  done
  Ncolumns=$((2 + 2*Nrates + Ncounter + Nstatus))
  echo >> ${outfile}

  #2. loop runs 

  for runnumber in ${runs[@]} ; do 
    filejob="${DIR}/cpass${pass}.job*.run${runnumber}.done"
    filemerge="${DIR}/merge.cpass${pass}.run${runnumber}.done"
    fileOCDB=$(grep /ocdb.log ${filemerge} | awk '{print $1}')
    if ! $(/bin/ls ${filemerge} &>/dev/null) ; then
      echo "${filemerge} does not exist!" >>${errfile}
      continue
    elif ! $(/bin/ls ${filejob} &>/dev/null) ; then
      echo "${filejob} does not exist!" >>${errfile}
      echo -n ${runnumber} >> ${outfile}
      for i in $(seq ${Ncolumns}) ; do 
        echo -n "-1" >> ${outfile}
      done
      echo >> ${outfile}
      continue
    fi
    echo -n ${runnumber} >> ${outfile}
    #pass0status= grep '/ocdb.log' ${filemerge} |  cut -d' ' -f2 | tr OK x1 | tr BAD xx0 | tr -d 'x'
    passStatus=$(grep '/ocdb.log' ${filemerge} | grep OK | wc -l)
    echo -n " ${passStatus}" >> ${outfile}
    qaStatus=$(grep '/mergeMakeOCDB.log' ${filemerge} | grep OK | wc -l)
    echo -n " ${qaStatus}" >> ${outfile}

    #fill OK/BAD rates
    for i in $(seq 0 $((${Nrates}-1))) ; do
      var1=$(grep "/${ratesString[${i}]}.log" ${filejob} | grep OK | wc -l)
      var2=$(grep "/${ratesString[${i}]}.log" ${filejob} | grep BAD | wc -l)

      if [[ ${ratesString[${i}]} == "stderr" ]] ; then
        var1=$(grep "stderr" ${filejob} | grep OK | wc -l)
        var2=$(grep "stderr" ${filejob} | grep "rec.log" | grep BAD | wc -l)
      fi
      echo -n " ${var1}" >> ${outfile}
      echo -n " ${var2}" >> ${outfile}
    done

    if [[ -f ${fileOCDB} ]] ; then
      #fill counter
      for i in $(seq 0 $((${Ncounter}-1))) ; do
        var1=$(grep Monalisa ${fileOCDB} | grep ${counterString[${i}]} | cut -f2)
        echo -n " ${var1:-"-1"}" >> ${outfile}
      done

      #fill status
      for i in $(seq 0 $((${Nstatus}-1))) ; do
        var1=$(grep "calibration status=" ${fileOCDB} | grep ${statusString[${i}]/Status/} | cut -d'=' -f2)
        echo -n " ${var1:-"-1"}" >> ${outfile}
      done
    fi
    echo >> ${outfile}
  done
  alilog_info  "[END] goMakeSummaryTree() with following parameters $*" 
  return 0
)

aliroot()
{
  args=("$@")
  if [[ -n ${useProfilingCommand} ]]; then
    profilerLogFile="cpu.txt"
    [[ "${args[@]}" =~ rec ]] && profilerLogFile="cpu_rec.txt"
    [[ "${args[@]}" =~ Calib ]] && profilerLogFile="cpu_calib.txt"
    echo running "${useProfilingCommand} aliroot ${args[@]} &> ${profilerLogFile}"
    ${useProfilingCommand} aliroot "${args[@]}" &> ${profilerLogFile}
  else
    #to prevent an infinite recursion use "command aliroot" to disable
    #aliases and functions
    echo running command aliroot "${args[@]}"
    command aliroot "${args[@]}"
  fi
  return 0
}

main "$@"
# vi:syntax=zsh
