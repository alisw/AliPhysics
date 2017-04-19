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
lastJobID=

# Those files are expected to be under $ALICE_PHYSICS/PWGPP/scripts. If a version in the current
# directory exists, it takes precedence.
sourceUtilities=( alilog4bash.sh utilities.sh )

main()
{
  #run in proper mode depending on the selection
  for scr in "${sourceUtilities[@]}"; do
    if [[ -e $scr ]]; then
      echo "Sourcing $scr from current directory" 1>&2
      source $scr false
    else
      echo "Sourcing $scr from AliPhysics" 1>&2
      source $ALICE_PHYSICS/PWGPP/scripts/$scr false
    fi
  done
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
    "CPass2") goCPass2 "$@";;
    "MakeLocalOCDBaccessConfig") goMakeLocalOCDBaccessConfig "$@";;
    "MergeCPass0") goMergeCPass0 "$@";;
    "MergeCPass1") goMergeCPass1 "$@";;
    "MergeCPass2") goMergeCPass2 "$@";;
    "MakeFilteredTrees") goMakeFilteredTrees "$@";;
    "MakeSummary") goMakeSummary "$@"
                   exit $? ;;
    "run") goSubmitMakeflow "$@"
           exit $? ;;
    "submit") goSubmitBatch "$@";;
    "test") goTest "$@";;
    "GenerateMakeflow") goGenerateMakeflow "$@";;
    "PrintValues") goPrintValues "$@";;
    "CreateQAplots") goCreateQAplots "$@";;
    "WaitForOutput") goWaitForOutput "$@";;
    "Merge") goMerge "$@";;
    "ppbench") goppbench "$@";;
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
  if [ -n "${SLURM_ARRAY_TASK_ID}"]; then
    SEED=${SLURM_ARRAY_JOB_ID}${SLURM_ARRAY_TASK_ID}
  fi
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

goppbench() (
  alilog_info "[BEGIN] goppbench() with the following extra parameters $*"
  configFile=$1
  shift
  parseConfig "configFile=$configFile" "$@" || return 1
  [[ "$SKIP_PPBENCH" == 1 ]] && { alilog_info "[END] goppbench() skipping ppbench"; \
                                  touch ppbench.done; return 0; }
  rm -rf AliRoot/ ppbench.done
  # Check if we have run it already. Prevent useless retries.
  maxCopyTries=2 xCopy -f -c -d . $commonOutputPath/ppbench/full_output.log
  [[ -e full_output.log ]] && { alilog_info "[END] goppbench() already ran, skipping"; return 0; }
  git clone http://git.cern.ch/pub/AliRoot --depth=1 || return 1
  pushd AliRoot
    export OCDB_TEST_ROOT=$PWD/OCDB
    test/runTests ppbench --debug --exit-on-error --variants default 2>&1 | tee -a full_output.log
    RV=${PIPESTATUS[0]}  # if this is nonzero we have a problem (caught later)
    pushd test/ppbench
      summarizeLogs * > summary_ppbench.log
      mv -v ../../full_output.log .
      xCopy -d $commonOutputPath/ppbench .
    popd
  popd
  [[ $RV != 0 ]] || cp -v AliRoot/test/ppbench/summary_ppbench.log ppbench.done
  rm -rf AliRoot/
  # Final output: <commonOutputPath>/ppbench.done                -- only upon success
  #               <commonOutputPath>/ppbench/full_output.log     -- raw full output
  #               <commonOutputPath>/ppbench/summary_ppbench.log -- same as ppbench.done on success
  # Validation will not continue (due to the lack of ppbench.done) in case of error
  alilog_info "[END] goppbench() finished"
)

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

goCPass2() (
  # Wrapper function that calls goCPass with the CPass1 option.
  declare -a params
  params=("$1" "$2" "$3" "$4" "$5" "$6" "$7" CPass2)
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
  cpass=$8  # cpass=CPass0 or CPass1 or CPass2
  shift 8
  extraOpts=("$@")
  cpass=${cpass##*CPass}

  parseConfig configFile=$configFile "$@" || return 1

  bigEcho "CPass ${cpass}"
  echo Start: goCPass${cpass}
  alilog_info "[BEGIN] goCPass${cpass}() with following extra parameters $*"

  # Remember working directory provided by the batch system (i.e. current dir).
  batchWorkingDirectory=$PWD

  # Use the jobindex only if set and non-negative
  if [[ -z "$jobindex" || "$jobindex" -lt 0 ]]; then
    [[ -n "$LSB_JOBINDEX" ]] && jobindex=$LSB_JOBINDEX
    [[ -n "$SGE_TASK_ID" ]] && jobindex=$SGE_TASK_ID
    [[ -n "$SLURM_ARRAY_TASK_ID" ]] && jobindex=$SLURM_ARRAY_TASK_ID
    # FIXME: should here also be a check for <0?

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

  logOutputDir=$runpath
  [[ -n "$logToFinalDestination" ]] && logOutputDir=${outputDir}
  if [[ -z "$dontRedirectStdOutToLog" ]]; then
    # Redirect all output to both file and console. Save fd 3 to restore later.
    exec 3>&1
    exec &> >(tee ${logOutputDir}/stdout)
  fi
  echo "$0 $*"

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
  echo "Making $infile visible in the current directory, $PWD"
  [[ $copyInputData == 0 ]] && ln -s ${infile} ${runpath}/${chunkName} \
                            || { maxCopyTries=1 remoteCpTimeout=1800 copyFileFromRemote ${infile} ${runpath}/ || exit 1; }

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
    0) filesCPassCustom=( "${batchWorkingDirectory}/runCPass0.sh"
                          "${batchWorkingDirectory}/recCPass0.C"
                          "${batchWorkingDirectory}/runCalibTrain.C"
                          "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                          "${batchWorkingDirectory}/localOCDB.tgz"
                          "${batchWorkingDirectory}/OCDB.root" )
       filesCPass=( "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/runCPass0.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/recCPass0.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/runCalibTrain.C" ) ;;

    1) filesCPassCustom=( "${batchWorkingDirectory}/runCPass1.sh"
                          "${batchWorkingDirectory}/recCPass1.C"
                          "${batchWorkingDirectory}/recCPass1_OuterDet.C"
                          "${batchWorkingDirectory}/runCalibTrain.C"
                          "${batchWorkingDirectory}/QAtrain_duo.C"
                          "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                          "${batchWorkingDirectory}/${configFile}"
                          "${batchWorkingDirectory}/OCDB.root" )
       filesCPass=( "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/runCalibTrain.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeQAgroups.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/runCPass1.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/recCPass1.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/recCPass1_OuterDet.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/QAtrain_duo.C"
                    "${commonOutputPath}/meta/cpass0.localOCDB.${runNumber}.tgz" ) ;;

    2) filesCPassCustom=( "${batchWorkingDirectory}/OCDB.root"
                          "${batchWorkingDirectory}/AODtrain.C"
                          "${batchWorkingDirectory}/rec.C"
                          "${batchWorkingDirectory}/raw2clust.C"
                          "${batchWorkingDirectory}/runPPass_pp.sh"
                          "${batchWorkingDirectory}/runPPass_pbpb.sh" )
       filesCPass=( "$ALICE_ROOT/test/QA/tag.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/PPass/AODtrain.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/PPass/rec.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/PPass/raw2clust.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/PPass/runPPass_pp.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/PPass/runPPass_pbpb.sh"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/QAtrain_duo.C"
                    "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeQAgroups.C"
                    "${commonOutputPath}/meta/cpass1.localOCDB.${runNumber}.tgz" ) ;;
  esac

  #first check if we have any custom scripts
  # -c: check if local source exists; -C: do not copy if local dest exists already
  # -f: copy all in the same dest dir (flat copy)
  xCopy -f -c -C -d . "${filesCPassCustom[@]}"
  for file in ${filesCPassCustom[*]}; do
    [[ ${file##*/} =~ .*\.sh ]] && printExec chmod +x ${file##*/}
  done

  #then download any missing ones from the default location
  # -c: check if local source exists; -C: do not copy if local dest exists already
  # -f: copy all in the same dest dir (flat copy)
  xCopy -f -c -C -d . "${filesCPass[@]}"
  for file in ${filesCPass[*]}; do
    [[ ${file##*/} =~ .*\.sh ]] && printExec chmod +x ${file##*/}
  done

  listDir "$PWD" "before running CPass${cpass}"

  # Monkey patching: emove spaces from around arguments to root macros. For example this sometimes
  # is known to fail: root 'macro.C(argument1, argument2)'
  sed -i '/.*root .*\.C/ s|\s*,\s*|,|g' *.sh

  # If OCDB is found here, then create a macro that configures local OCDB access.
  # This step also decompresses the tarball into $PWD/OCDB.

  # this would be relevant only for cpass0 :
  # custom initial specific OCDB objects provided by the user at the beginning
  [[ -f localOCDB.tgz && $cpass == 0 ]] && tar xzvvf localOCDB.tgz

  ocdbTarball=cpass$(($cpass-1)).localOCDB.${runNumber}.tgz
  [[ -f $ocdbTarball ]] \
    && printExec goMakeLocalOCDBaccessConfig $ocdbTarball \
    || echo "WARNING: file $ocdbTarball not found!"

  # A post-setup action (different for CPass0/1) can be defined. Here it is executed.
  case $cpass in
    0) postSetUpActionCPass="$postSetUpActionCPass0" ;;
    1) postSetUpActionCPass="$postSetUpActionCPass1" ;;
    2) postSetUpActionCPass="$postSetUpActionCPass2" ;;
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
        for fakeOutput in AliESDs.root AliESDfriends.root CalibObjects.root rec.log calib.log \
            syswatch.log syswatch_rec.log; do
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

      # Check if previous pass has produced calibration.
      if [[ ! $(/bin/ls -1 OCDB/*/*/*/*.root 2>/dev/null) ]]; then
        touch $doneFileTmp
        echo "CPass$(($cpass-1)) produced no calibration! Exiting..." >> $doneFileTmp
        copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
        return 1
      fi

      # Create the Barrel and OuterDet directories for CPass1 and link the local OCDB directory
      # there to make the localOCDBaccessConfig.C file work, since it may point to the OCDB
      # entries using a relative path, e.g. local://./OCDB.
      echo "Linking the OCDB/ for Barrel and OuterDet directories"
      mkdir Barrel OuterDet
      ln -s ../OCDB Barrel/OCDB
      ln -s ../OCDB OuterDet/OCDB
      listDir "$PWD" "after linking the OCDB for Barrel and OuterDet"

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
        for fakeOutput in AliESDs_Barrel.root AliESDfriends_Barrel.root CalibObjects.root \
                          QAresults_barrel.root EventStat_temp_barrel.root AODtpITS.root \
                          AliESDs_Outer.root AliESDfriends_Outer.root QAresults_outer.root \
                          EventStat_temp_outer.root rec.log calib.log qa.log filtering.log \
                          syswatch_rec_Barrel.log syswatch_rec_Outer.log syswatch_calib.log \
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
      [[ ! -f CalibObjects.root && -f Barrel/CalibObjects.root ]] && mv Barrel/CalibObjects.root .
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
    2) 
      # Check if previous pass has produced calibration.
      if [[ ! $(/bin/ls -1 OCDB/*/*/*/*.root 2>/dev/null) ]] && [ $runCPass1MergeMakeOCDB -eq 1 ] ; then
        touch $doneFileTmp
        echo "CPass$(($cpass-1)) produced no calibration! Exiting..." >> $doneFileTmp
        copyFileToRemote "$doneFileTmp" "$(dirname "$doneFile")" || rm -f "$doneFileTmp"
        return 1
      fi

      if [[ -n "$pretend" ]]; then
        echo "Pretending to run: $runpath/runCPass0.sh /$infile $nEvents $runNumber $ocdbPath $recoTriggerOptions"
        sleep $pretendDelay
        for fakeOutput in AliESDs.root AliESDfriends.root \
                          QAresults.root QAresults.root \
                          EventStat_temp_outer.root rec.log qa.log filtering.log \
                          syswatch_rec.log AliAOD.root; do
          touch $fakeOutput
        done
      else
        collisionSystem=$(run2collisionSystem "$runNumber")
        printExec ./runPPass_${collisionSystem}.sh "/$infile" "SPLIT" "$nEvents" "$runNumber" "$ocdbPath"
	[[ -f AliESDs.root ]] && echo AliESDs.root > filtered.list
	printExec goMakeFilteredTrees $PWD $runNumber "$PWD/filtered.list" $filteringFactorHighPt \
                            $filteringFactorV0s $ocdbPath 1000000 0 10000000 0 \
                            $configFile AliESDs.root "${extraOpts[@]}"
      fi
    ;;

  esac

  listDir "$PWD" "after running CPass${cpass}"

  # CPass has completed. Copy all created files to the destination (which might be remote).
  # Note: stdout is not copied, this will happen at the very end.
  printExec rm -f ./$chunkName
  filesToCopy=()
  while read cpdir; do
    filesToCopy+=($cpdir/!(stdout|cpass$(($cpass-1))*.tgz))
  done < <(find . -type d)
  xCopy -d $outputDir/ "${filesToCopy[@]}"

  # Validate CPass.
  case $cpass in

    0)
      # Validate CPass0.
      if summarizeLogs * */* | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
        [[ -f $outputDirMC/galice.root ]] && echo "sim $outputDirMC/galice.root" >> $doneFileTmp
        [[ -f AliESDfriends_v1.root ]] && echo "calibfile $outputDir/AliESDfriends_v1.root" >> $doneFileTmp
        [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" >> $doneFileTmp  # new name of AliESDfriends_v1.root
        [[ -f AliESDs.root ]] && echo "esd $outputDir/AliESDs.root" >> $doneFileTmp
      fi
      #report syswatch logs
      reportDoneFile syswatchRec syswatch_rec.log ${outputDir} >> $doneFileTmp
      reportDoneFile syswatchCalib syswatch.log ${outputDir} >> $doneFileTmp

      # End of CPass0 validation.
    ;;

    1)
      # Validate CPass1.
      if summarizeLogs * */* | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
        [[ -f AliESDs_Barrel.root ]] && echo "esd $outputDir/AliESDs_Barrel.root" >> $doneFileTmp
        [[ -f AliESDfriends_v1.root ]] && echo "calibfile $outputDir/AliESDfriends_v1.root" >> $doneFileTmp
        [[ -f CalibObjects.root ]] && echo "calibfile $outputDir/CalibObjects.root" >> $doneFileTmp  # new name of AliESDfriends_v1.root
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
      #report syswatch logs
      reportDoneFile syswatchRec syswatch_rec_Outer.log ${outputDir} >> $doneFileTmp
      reportDoneFile syswatchRec syswatch_rec_Barrel.log ${outputDir} >> $doneFileTmp
      reportDoneFile syswatchCalib syswatch_calib.log ${outputDir} >> $doneFileTmp
      # End of CPass1 validation.
    ;;
    2)
      #validate CPass2
      if summarizeLogs * */* | sed -e "s|$PWD|$outputDir|" >> $doneFileTmp; then
        reportDoneFile esd AliESDs.root $outputDir >> $doneFileTmp
        reportDoneFile aod AliAOD.root $outputDir >> $doneFileTmp
        reportDoneFile syswatchRec syswatch_rec.log ${outputDir} >> $doneFileTmp
        reportDoneFile qafile QA_results.root $outputDir >> $doneFileTmp
        reportDoneFile qafile QAresults.root $outputDir >> $doneFileTmp
      fi
    ;;

  esac
  echo "dir $outputDir" >> $doneFileTmp
  
  # Copy stdout to destination.
  if [[ -z "$dontRedirectStdOutToLog" ]]; then
    echo "Copying stdout. NOTE: this is the last bit you will see in the log!"
    exec 1>&3 3>&-
    exec 2>&1
    copyFileToRemote stdout $outputDir
  fi

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

goMergeCPass2() (
  goMergeCPass CPass2 "$@"
)

goMergeCPass()
(
  # Find CPass0 output files and merge them. OCDB is created.

  # First argument can be either CPass0 or CPass1.
  cpass=$1
  cpass=${cpass##*CPass}
  [[ $cpass != 0 && $cpass != 1 && $cpass != 2 ]] \
    && echo "FATAL: call goMergeCPass with CPass[0|1] as first param!" && return 1
  shift

  # Arguments in common between MergeCPass0 and MergeCPass1.
  outputDir=$1
  ocdbStorage=$2
  configFile=$3
  export runNumber=$4
  shift 4

  #some defaults to avoid empty strings
  calibrationFilesToMerge="dummy"
  qaFilesToMerge="dummy"
  filteredFilesToMerge="dummy"
  syslogsRecToMerge="dummy"
  syslogsCalibToMerge="dummy"

  parseConfig configFile=$configFile "$@" || return 1

  bigEcho "Merging CPass ${cpass}"
  echo Start: goMergeCPass${cpass}
  alilog_info "[BEGIN] goMergeCPass${cpass}() with following parameters $*"

  listDir "$PWD" "before transferring files locally in MergeCPass${cpass}"

  # Copy all the files to a local dir tree. Replace remote file names in the lists with their local
  # versions.
  for remoteList in $calibrationFilesToMerge \
                    $qaFilesToMerge \
                    $filteredFilesToMerge \
                    $syslogsRecToMerge \
                    $syslogsCalibToMerge; do
    localList=local.${remoteList##*/}
    rm -f "$localList" && touch "$localList"
    while read sourceFile; do
      destinationFile="${PWD}/${sourceFile#${commonOutputPath}}"
      copyFileFromRemote "$sourceFile" "$(dirname "${destinationFile}")" && \
        echo "$destinationFile" >> "$localList"
    done < <(cat "$remoteList")
  done
  # Subsequent calls will use the local.* version of those files.
  calibrationFilesToMerge=local.${calibrationFilesToMerge##*/}
  filteredFilesToMerge=local.${filteredFilesToMerge##*/}
  syslogsRecToMerge=local.${syslogsRecToMerge##*/}
  syslogsCalibToMerge=local.${syslogsCalibToMerge##*/}

  if [[ "$qaFilesToMerge" != '' ]]; then
    qaFilesToMerge=local.${qaFilesToMerge##*/}
    #strip filenames as QA merging requires list of directories
    # Important to have the string "Stage.txt" in the filename to trigger the merging.
    # It has to be a list of directories containing the files.
    sed -e 's|/[^\/]*$||g' "$qaFilesToMerge" > ${qaFilesToMerge}.lastMergingStage.txt.list
    qaFilesToMerge=${qaFilesToMerge}.lastMergingStage.txt.list
  fi

  # Record the working directory provided by the batch system.
  batchWorkingDirectory=$PWD

  [[ -z "$commonOutputPath" ]] && commonOutputPath=$PWD

  # This file signals that the job is done, no matter the result
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
  if [[ -z "$dontRedirectStdOutToLog" ]]; then
    # Redirect all output to both file and console. Save fd 3 to restore later.
    exec 3>&1
    exec &> >(tee ${logOutputDir}/stdout)
  fi
  echo "$0 $*"

  mergingScript="mergeMakeOCDB.byComponent.perStage.sh"

  if [[ $cpass -ge 1 ]]; then
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
  
  # Copy files in case they are not already there.
  case $cpass in
    0) filesMergeCPassCustom=( "${batchWorkingDirectory}/${calibrationFilesToMerge}"
                               "${batchWorkingDirectory}/${syslogsRecToMerge}"
                               "${batchWorkingDirectory}/${syslogsCalibToMerge}"
                               "${batchWorkingDirectory}/${mergingScript}"
                               "${batchWorkingDirectory}/OCDB.root"
                               "${batchWorkingDirectory}/localOCDB.tgz"
                               "${batchWorkingDirectory}/localOCDBaccessConfig.C" )
       filesMergeCPass=( "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/${mergingScript}"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/mergeByComponent.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/makeOCDB.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass0/merge.C" ) ;;

    1) filesMergeCPassCustom=( "${batchWorkingDirectory}/${calibrationFilesToMerge}"
                               "${batchWorkingDirectory}/${qaFilesToMerge}"
                               "${batchWorkingDirectory}/${filteredFilesToMerge}"
                               "${batchWorkingDirectory}/${syslogsRecToMerge}"
                               "${batchWorkingDirectory}/${syslogsCalibToMerge}"
                               "${batchWorkingDirectory}/${mergingScript}"
                               "${batchWorkingDirectory}/OCDB.root"
                               "${batchWorkingDirectory}/localOCDBaccessConfig.C" )
      filesMergeCPass=( "${commonOutputPath}/meta/cpass0.localOCDB.${runNumber}.tgz"
                         "${batchWorkingDirectory}/QAtrain_duo.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/${mergingScript}"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeByComponent.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/makeOCDB.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/merge.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/QAtrain_duo.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeQAgroups.C" ) ;;

    2) filesMergeCPassCustom=( "${batchWorkingDirectory}/${calibrationFilesToMerge}"
                               "${batchWorkingDirectory}/${qaFilesToMerge}"
                               "${batchWorkingDirectory}/${filteredFilesToMerge}"
                               "${batchWorkingDirectory}/${syslogsRecToMerge}"
                               "${batchWorkingDirectory}/${syslogsCalibToMerge}"
                               "${batchWorkingDirectory}/OCDB.root"
                               "${batchWorkingDirectory}/localOCDBaccessConfig.C"
                               "${batchWorkingDirectory}/QAtrain_duo.C" )
       filesMergeCPass=( "${commonOutputPath}/meta/cpass1.localOCDB.${runNumber}.tgz"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/mergeQAgroups.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/merge.C"
                         "${ALICE_PHYSICS}/PWGPP/CalibMacros/CPass1/QAtrain_duo.C" ) ;;
  esac

  #first check if we have any custom scripts
  # -c: check if local source exists; -C: do not copy if local dest exists already
  # -f: copy all in the same dest dir (flat copy)
  xCopy -f -c -C -d . "${filesMergeCPassCustom[@]}"
  for file in ${filesMergeCPassCustom[*]}; do
    [[ ${file##*/} =~ .*\.sh ]] && printExec chmod +x ${file##*/}
  done

  #then download any missing ones from the default location
  # -c: check if local source exists; -C: do not copy if local dest exists already
  # -f: copy all in the same dest dir (flat copy)
  xCopy -f -c -C -d . "${filesMergeCPass[@]}"
  for file in ${filesMergeCPass[*]}; do
    [[ ${file##*/} =~ .*\.sh ]] && printExec chmod +x ${file##*/}
  done

  # Monkey patching: remove spaces from around arguments to root macros. For example this sometimes
  # is known to fail: root 'macro.C(argument1, argument2)'
  sed -i '/.*root .*\.C/ s|\s*,\s*|,|g' *.sh

  alirootInfo > ALICE_ROOT.log

  # this would be relevant only for cpass0 :
  # custom initial specific OCDB objects provided by the user at the beginning
  [[ -f localOCDB.tgz && $cpass == 0 ]] && tar xzvvf localOCDB.tgz

  # Configure local OCDB storage (a macro is produced). Only used in MergeCPass1 but could be used
  # in CPass0 too. It's harmless in any case.
  ocdbTarball=cpass$(($cpass-1)).localOCDB.${runNumber}.tgz
  if [[ -f $ocdbTarball ]]; then
    printExec goMakeLocalOCDBaccessConfig "$ocdbTarball"
  elif [[ $cpass -ge 1 ]]; then
    # Print a warning only on CPass1.
    echo "WARNING: file $ocdbTarball not found!"
  fi

  listDir ".." "before running MergeCPass${cpass} in ${PWD}"

  # Merge calibration.
  chmod u+x $mergingScript
  mkdir -p OCDB

  # Run Merge of this CPass.
  case $cpass in
    0)
      # MergeCPass0.

      if [[ -n ${pretend} ]]; then
        sleep $pretendDelay
        for file in CalibObjects.root ocdb.log merge.log dcsTime.root syswatch.log; do
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
    2)
      if [[ -n ${pretend} ]]; then
        sleep ${pretendDelay}
        for file in "${qaMergedOutputFileName}" \
                    mergeQA.log syswatch.log; do
          touch $file
        done
      else
        printExec aliroot -b -q "QAtrain_duo.C(\"\",${runNumber},\"${qaFilesToMerge}\",1,\"${ocdbStorage}\")" > mergeQA.log
        #QAtrain_duo default output is called QAresults.root, move to the expected name
        mv QAresults.root ${qaMergedOutputFileName}
      fi
    ;;
  esac
  
  # Create tarball with OCDB, store on the shared directory, create signal file on batch directory
  baseTar="cpass${cpass}.localOCDB.${runNumber}.tgz"
  printExec tar czvvf ${batchWorkingDirectory}/${baseTar} ./OCDB && \
    echo "ocdbTarball $commonOutputPath/meta/$baseTar" >> $doneFileTmp

  listDir "$batchWorkingDirectory" "after tarball creation in MergeCPass${cpass}"

  #merge the syslogs
  echo "trying to merge syslogs $syslogsRecToMerge $syslogsCalibToMerge"
  echo "in dir $PWD"
  if [[ -f "$syslogsRecToMerge" ]]; then
    echo "mergeSysLogs syswatch.rec.cpass${cpass}.tree @${syslogsRecToMerge}"
    mergeSysLogs syswatch.rec.cpass${cpass}.tree @${syslogsRecToMerge}
  fi
  if [[ -f "$syslogsCalibToMerge" ]]; then
    echo "mergeSysLogs syswatch.calib.cpass${cpass}.tree @${syslogsCalibToMerge}"
    mergeSysLogs syswatch.calib.cpass${cpass}.tree @${syslogsCalibToMerge}
  fi

  # CPass has completed. Copy all created files to the destination (which might be remote).
  # Note: stdout is not copied, this will happen at the very end.
  filesToCopy=()
  while read cpdir; do
    filesToCopy+=($cpdir/!(stdout))
  done < <(find . -type d)
  xCopy -d $outputDir/ "${filesToCopy[@]}"

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
  reportDoneFile syswatchRec syswatch.rec.cpass${cpass}.tree ${outputDir} >> $doneFileTmp
  reportDoneFile syswatchCalib syswatch.calib.cpass${cpass}.tree ${outputDir} >> $doneFileTmp

  # Copy stdout to destination.
  if [[ -z "$dontRedirectStdOutToLog" ]]; then
    echo "Copying stdout. NOTE: this is the last bit you will see in the log!"
    exec 1>&3 3>&-
    exec 2>&1
    copyFileToRemote stdout $outputDir
  fi

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

  # Copy to the current directory the extra library scripts, if not there.
  for scr in "${sourceUtilities[@]}"; do
    [[ -e $scr ]] || printExec cp -v $ALICE_PHYSICS/PWGPP/scripts/$scr .
  done

  self=$0

  # For reference copy the setup to the output dir.
  copyFileToRemote $self $commonOutputPath
  copyFileToRemote $configFile $commonOutputPath
  copyFileToRemote $inputFileList $commonOutputPath
  copyFileToRemote "${sourceUtilities[@]}" $commonOutputPath

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
  xCopy -f -C -d $PWD $commonOutputPath/benchmark.makeflow.makeflowlog  # do not copy if it exists locally
  awk '/STARTED/   {startTime=$3}
       /COMPLETED/ {endTime=$3}
       END         {print "makeflow running time: "(endTime-startTime)/1000000/3600" hours"}' \
      benchmark.makeflow.makeflowlog > running_time
  cat running_time >> summary.log
  cat running_time >> summary_full.log
  xCopy -f -d $commonOutputPath summary.log summary_full.log running_time

  goCheckSummary
  return $?
}

goCheckSummary() {
  # Checks if summary.log (which must be available in the current dir) has any
  # indication of errors.
  local BAD
  [[ ! -f summary.log ]] && return 1
  BAD=$(awk '/error summary:/,/detailed summary:/' summary.log | grep -v '===' | sed -e '/^$/d' | wc -l)
  BAD=$((BAD && 1))
  [[ $BAD == 0 ]] && printf "\n\n==> No BAD files found!\n\n\n" \
                  || printf "\n\n==> Validation failed: some BAD files were found! Check summary.log for more details.\n\n\n"
  return $BAD
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
  inputFiles=( OCDB.root
               localOCDB.tgz
               localOCDBaccessConfig.C
               QAtrain_duo.C
               runCPass1.sh
               recCPass1.C
               recCPass1_OuterDet.C
               runCalibTrain.C
               runCPass0.sh
               recCPass0.C
               runQA.sh
               AODtrain.C
               rec.C
               raw2clust.C
               runPPass_pp.sh
               runPPass_pbpb.sh
               AODtrain.C
               mergeQAgroups.C )
  for file in ${inputFiles[*]}; do
    [[ -f ${file} ]] && copyFiles+=("${file}")
  done

  #create the makeflow file
  [[ -n ${batchFlags} ]] && echo "BATCH_OPTIONS = ${batchFlags}"
  declare -A arr_cpass0_merged arr_cpass1_merged
  declare -A arr_cpass0_calib_list arr_cpass1_calib_list 
  declare -A arr_cpass1_QA_files_list arr_cpass1_ESD_list arr_cpass1_filtered_list
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

    #ppbench
    echo "### ppbench ###"
    echo "ppbench.done: benchmark.sh ${sourceUtilities[*]} ${configFile}"
    echo -e "\t${alirootEnv} ./benchmark.sh ppbench ${configFile} ${extraOpts[@]}"" "
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
      echo "${arr_cpass0_outputs[${jobindex}]}: benchmark.sh ppbench.done ${sourceUtilities[*]} ${configFile} ${copyFiles[@]}"
      echo -e "\t${alirootEnv} ./benchmark.sh CPass0 \$OUTPATH/000${runNumber}/cpass0 ${inputFile} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]}"" "
      echo ; echo

      #CPass1
      #arr_cpass1_outputs[${jobindex}]="${commonOutputPath}/meta/cpass1.job${jobindex}.run${runNumber}.done"
      arr_cpass1_outputs[${jobindex}]="cpass1.job${jobindex}.run${runNumber}.done"
      echo "### CPass1 ###"
      echo "${arr_cpass1_outputs[${jobindex}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} merge.cpass0.run${runNumber}.done ${copyFiles[@]}"
      echo -e "\t${alirootEnv} ./benchmark.sh CPass1 \$OUTPATH/000${runNumber}/cpass1 ${inputFile} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]}"" "
      echo ; echo

      #CPass2
      #arr_cpass2_outputs[${jobindex}]="${commonOutputPath}/meta/cpass2.job${jobindex}.run${runNumber}.done"
      arr_cpass2_outputs[${jobindex}]="cpass2.job${jobindex}.run${runNumber}.done"
      echo "### CPass2 ###"
      echo "${arr_cpass2_outputs[${jobindex}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} merge.cpass1.run${runNumber}.done ${copyFiles[@]}"
      echo -e "\t${alirootEnv} ./benchmark.sh CPass2 \$OUTPATH/000${runNumber}/cpass2 ${inputFile} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]}"" "
      echo ; echo
      ((jobindex++))

    done< <(grep "/000${runNumber}/" ${inputFileList})

    #######################CPass0############################

    #CPass0 list of Calib files to merge
    #arr_cpass0_calib_list[${runNumber}]="${commonOutputPath}/meta/cpass0.calib.run${runNumber}.list"
    arr_cpass0_calib_list[${runNumber}]="cpass0.calib.run${runNumber}.list"
    echo "### Produces the list of CPass0 files to merge (executes locally) ###"
    echo "${arr_cpass0_calib_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass0_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues calibfile ${arr_cpass0_calib_list[${runNumber}]} ${arr_cpass0_outputs[*]}"
    echo ; echo

    #CPass0 list of rec syslogs to merge
    arr_cpass0_rec_syswatch_list[${runNumber}]="cpass0.syswatch.rec.run${runNumber}.list"
    echo "### Produces the list of CPass0 files to merge (executes locally) ###"
    echo "${arr_cpass0_rec_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass0_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchRec ${arr_cpass0_rec_syswatch_list[${runNumber}]} ${arr_cpass0_outputs[*]}"
    echo ; echo

    #CPass0 list of calib syslogs to merge
    arr_cpass0_calib_syswatch_list[${runNumber}]="cpass0.syswatch.calib.run${runNumber}.list"
    echo "### Produces the list of CPass0 files to merge (executes locally) ###"
    echo "${arr_cpass0_calib_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass0_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchCalib ${arr_cpass0_calib_syswatch_list[${runNumber}]} ${arr_cpass0_outputs[*]}"
    echo ; echo

    #CPass0 merging
    echo "### Merges CPass0 files ###"
    #arr_cpass0_merged[${runNumber}]="${commonOutputPath}/meta/merge.cpass0.run${runNumber}.done"
    arr_cpass0_merged[${runNumber}]="merge.cpass0.run${runNumber}.done"
    echo "${arr_cpass0_merged[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} ${arr_cpass0_calib_list[${runNumber}]} ${arr_cpass0_rec_syswatch_list[${runNumber}]} ${arr_cpass0_calib_syswatch_list[${runNumber}]} ${copyFiles[@]}"
    echo -e "\t${alirootEnv} ./benchmark.sh MergeCPass0 \$OUTPATH/000${runNumber}/cpass0 ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${arr_cpass0_calib_list[${runNumber}]} syslogsRecToMerge=${arr_cpass0_rec_syswatch_list[${runNumber}]} syslogsCalibToMerge=${arr_cpass0_calib_syswatch_list[${runNumber}]} ${extraOpts[@]}"" "
    echo ; echo

    #######################CPass1############################

    # CPass1 list of QA files.
    arr_cpass1_QA_files_list[${runNumber}]="cpass1.QA.run${runNumber}.list"
    echo "### Lists CPass1 QA files ###"
    echo "${arr_cpass1_QA_files_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues qafile ${arr_cpass1_QA_files_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    # CPass1 list of calib files
    #arr_cpass1_calib_list[${runNumber}]="${commonOutputPath}/meta/cpass1.calib.run${runNumber}.list"
    arr_cpass1_calib_list[${runNumber}]="cpass1.calib.run${runNumber}.list"
    echo "### Lists CPass1 Calib ###"
    echo "${arr_cpass1_calib_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues calibfile ${arr_cpass1_calib_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    # CPass1 list of filtered tree files
    #arr_cpass1_filtered_list[${runNumber}]="${commonOutputPath}/meta/cpass1.filtered.run${runNumber}.list"
    arr_cpass1_filtered_list[${runNumber}]="cpass1.filtered.run${runNumber}.list"
    echo "### Lists CPass1 filtered ###"
    echo "${arr_cpass1_filtered_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues filteredTree ${arr_cpass1_filtered_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    #CPass1 list of rec syslogs to merge
    arr_cpass1_rec_syswatch_list[${runNumber}]="cpass1.syswatch.rec.run${runNumber}.list"
    echo "### Produces the list of CPass1 files to merge (executes locally) ###"
    echo "${arr_cpass1_rec_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchRec ${arr_cpass1_rec_syswatch_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    #CPass1 list of calib syslogs to merge
    arr_cpass1_calib_syswatch_list[${runNumber}]="cpass1.syswatch.calib.run${runNumber}.list"
    echo "### Produces the list of CPass1 files to merge (executes locally) ###"
    echo "${arr_cpass1_calib_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass1_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchCalib ${arr_cpass1_calib_syswatch_list[${runNumber}]} ${arr_cpass1_outputs[*]}"
    echo ; echo

    #CPass1 merging
    #arr_cpass1_merged[${runNumber}]="${commonOutputPath}/meta/merge.cpass1.run${runNumber}.done"
    arr_cpass1_merged[${runNumber}]="merge.cpass1.run${runNumber}.done"
    echo "### Merges CPass1 files ###"
    echo "${arr_cpass1_merged[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} ${arr_cpass1_calib_list[${runNumber}]} ${arr_cpass1_QA_files_list[${runNumber}]} ${arr_cpass1_filtered_list[${runNumber}]} ${copyFiles[@]}"
    echo -e "\t${alirootEnv} ./benchmark.sh MergeCPass1 \$OUTPATH/000${runNumber}/cpass1 ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${arr_cpass1_calib_list[${runNumber}]} qaFilesToMerge=${arr_cpass1_QA_files_list[${runNumber}]} filteredFilesToMerge=${arr_cpass1_filtered_list[${runNumber}]} syslogsRecToMerge=${arr_cpass1_rec_syswatch_list[${runNumber}]} syslogsCalibToMerge=${arr_cpass1_calib_syswatch_list[${runNumber}]} ${extraOpts[@]}"
    echo ; echo

    #######################CPass2############################

    # CPass2 list of QA files.
    arr_cpass2_QA_files_list[${runNumber}]="cpass2.QA.run${runNumber}.list"
    echo "### Lists CPass2 QA files ###"
    echo "${arr_cpass2_QA_files_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass2_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues qafile ${arr_cpass2_QA_files_list[${runNumber}]} ${arr_cpass2_outputs[*]}"
    echo ; echo

    # CPass2 list of calib files
    #arr_cpass2_calib_list[${runNumber}]="${commonOutputPath}/meta/cpass2.calib.run${runNumber}.list"
    arr_cpass2_calib_list[${runNumber}]="cpass2.calib.run${runNumber}.list"
    echo "### Lists CPass2 Calib ###"
    echo "${arr_cpass2_calib_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass2_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues calibfile ${arr_cpass2_calib_list[${runNumber}]} ${arr_cpass2_outputs[*]}"
    echo ; echo

    # CPass2 list of filtered tree files
    #arr_cpass2_filtered_list[${runNumber}]="${commonOutputPath}/meta/cpass2.filtered.run${runNumber}.list"
    arr_cpass2_filtered_list[${runNumber}]="cpass2.filtered.run${runNumber}.list"
    echo "### Lists CPass2 filtered ###"
    echo "${arr_cpass2_filtered_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass2_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues filteredTree ${arr_cpass2_filtered_list[${runNumber}]} ${arr_cpass2_outputs[*]}"
    echo ; echo

    #CPass2 list of rec syslogs to merge
    arr_cpass2_rec_syswatch_list[${runNumber}]="cpass2.syswatch.rec.run${runNumber}.list"
    echo "### Produces the list of CPass2 files to merge (executes locally) ###"
    echo "${arr_cpass2_rec_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass2_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchRec ${arr_cpass2_rec_syswatch_list[${runNumber}]} ${arr_cpass2_outputs[*]}"
    echo ; echo

    #CPass2 list of calib syslogs to merge
    arr_cpass2_calib_syswatch_list[${runNumber}]="cpass2.syswatch.calib.run${runNumber}.list"
    echo "### Produces the list of CPass2 files to merge (executes locally) ###"
    echo "${arr_cpass2_calib_syswatch_list[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${arr_cpass2_outputs[*]}"
    echo -e "\tLOCAL ./benchmark.sh PrintValues syswatchCalib ${arr_cpass2_calib_syswatch_list[${runNumber}]} ${arr_cpass2_outputs[*]}"
    echo ; echo

    #CPass2 merging
    #arr_cpass2_merged[${runNumber}]="${commonOutputPath}/meta/merge.cpass2.run${runNumber}.done"
    arr_cpass2_merged[${runNumber}]="merge.cpass2.run${runNumber}.done"
    echo "### Merges CPass2 files ###"
    echo "${arr_cpass2_merged[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} ${arr_cpass2_calib_list[${runNumber}]} ${arr_cpass2_QA_files_list[${runNumber}]} ${arr_cpass2_filtered_list[${runNumber}]} ${arr_cpass2_calib_syswatch_list[${runNumber}]} ${arr_cpass2_rec_syswatch_list[${runNumber}]} ${copyFiles[@]}"
    echo -e "\t${alirootEnv} ./benchmark.sh MergeCPass2 \$OUTPATH/000${runNumber}/cpass2 ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${arr_cpass2_calib_list[${runNumber}]} qaFilesToMerge=${arr_cpass2_QA_files_list[${runNumber}]} filteredFilesToMerge=${arr_cpass2_filtered_list[${runNumber}]} syslogsRecToMerge=${arr_cpass2_rec_syswatch_list[${runNumber}]} syslogsCalibToMerge=${arr_cpass2_calib_syswatch_list[${runNumber}]} ${extraOpts[@]}"
    echo ; echo

  ############################# DEBUG/profiling ##################################

    #CPass0 wrapped in a profiling tool (valgrind,....)
    if [[ -n ${profilingCommand} ]]; then
      inputFile=$(grep -m1 "${runNumber}/" ${inputFileList})
      [[ -z ${nEventsProfiling} ]] && nEventsProfiling=2
      currentDefaultOCDB=$(setYear ${inputFile} ${defaultOCDB})
      jobindex="profiling"

      #arr_cpass0_profiled_outputs[${runNumber}]="${commonOutputPath}/meta/profiling.cpass0.job${jobindex}.run${runNumber}.done"
      arr_cpass0_profiled_outputs[${runNumber}]="profiling.cpass0.job${jobindex}.run${runNumber}.done"
      echo "### CPass0 in a profiler ###"
      echo "${arr_cpass0_profiled_outputs[${runNumber}]}: benchmark.sh ${sourceUtilities[*]} ${configFile} ${copyFiles[@]}"
      profilingCommand=$(encSpaces "${profilingCommand}")
      echo -e "\t${alirootEnv} ./benchmark.sh CPass0 \$OUTPATH/000${runNumber}/${jobindex} ${inputFile} ${nEventsProfiling} ${currentDefaultOCDB} ${configFile} ${runNumber} ${jobindex} ${extraOpts[@]} useProfilingCommand=${profilingCommand}"
      echo ; echo
    fi

  done #runs

  ############################# Summary ##################################

  echo "### Summary ###"
  echo "summary.log: benchmark.sh ${sourceUtilities[*]} ${configFile} ${arr_cpass0_outputs[*]} ${arr_cpass0_merged[*]} ${arr_cpass1_outputs[*]} ${arr_cpass1_merged[*]} ${arr_cpass2_outputs[*]} ${arr_cpass2_merged[*]}"
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

  [[ -e utilities.sh ]] && cp -f utilities.sh "$outputDir"

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
  echo   based on OCDB: ${localOCDBpathCPass0}
  echo

  local tempLocalOCDB=""
  if [[ -f localOCDBaccessConfig.C ]]; then
    tempLocalOCDB=$(mktemp -t tempLocalOCDB.XXXXXX)
    egrep "SetSpecificStorage" localOCDBaccessConfig.C > ${tempLocalOCDB}
  fi

  echo "localOCDBaccessConfig()"                               >  localOCDBaccessConfig.C
  echo "{"                                                     >> localOCDBaccessConfig.C
  echo "  AliCDBManager* man = AliCDBManager::Instance();"     >> localOCDBaccessConfig.C
  [[ -f "${tempLocalOCDB}" ]] && cat ${tempLocalOCDB}              >> localOCDBaccessConfig.C
  spitOutLocalOCDBaccessConfig ${localOCDBpathCPass0}|sort|uniq  >> localOCDBaccessConfig.C
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

  # steer batch farm type via batchCommand in config file
  batchSystem="SGE"
  if [[ "$batchCommand" =~ bsub ]]; then
    batchSystem="LSF"
  elif [[ "$batchCommand" =~ sbatch ]]; then
    batchSystem="SLURM"
  fi

  newFarm=$(which qsub|grep "^/usr/bin/qsub")
  
  #if [[ -z "${newFarm}" ]]
  if [ "$batchSystem" == "LSF" ]
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
  elif [ "$batchSystem" == "SGE" ]; then
    #new SGE farm
    if [[ ${waitForJOBID} =~ "000" ]]; then
      echo ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
      ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
    else
      echo ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -hold_jid "${waitForJOBID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
      ${batchCommand} ${batchFlags} -wd ${commonOutputPath} -b y -v commonOutputPath -N "${JobID}" -t "${startID}-${endID}" -hold_jid "${waitForJOBID}" -e "${commonOutputPath}/logs/" -o "${commonOutputPath}/logs/" "${command}" "${commandArgs[@]}"
    fi
  elif [ "$batchSystem" == "SLURM" ]; then
    # ---| in case a wait job is requested add it to the job arguments |-----------
    jobArgs=""
    if [[ ! ${waitForJOBID} =~ "000" ]]; then
      jobArgs="-d afterany:${waitForJOBID}"
    fi

    # ---| put together the job command and execute it, the job id will be saved |-
    jobcmd="${batchCommand} ${batchFlags} \
            -J ${JobID} --array=${startID}-${endID} \
            -o ${commonOutputPath}/logs/${JobID}.%j.%a.out -e ${commonOutputPath}/logs/${JobID}.%j.%a.err \
            --export=commonOutputPath \
            --workdir=${commonOutputPath} \
            ${jobArgs} \
            ${command} ${commandArgs[@]}"
    echo $jobcmd
    lastJobID=$(eval $jobcmd | tee /dev/tty | awk '{print $4}')
  else
    echo "unknown batch system: $batchSystem"
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
    JOBID6="p2_${JOBpostfix}"
    JOBID6wait="w2_${JOBpostfix}"
    JOBID7="m2_${JOBpostfix}"
    JOBID7wait="wm2_${JOBpostfix}"
    JOBID8="s1_${JOBpostfix}"
    JOBID8wait="ws1_${JOBpostfix}"
    JOBID9="QA_${JOBpostfix}"
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
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID1wait}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

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
                        "${configPath}/mergeMakeOCDB.byComponent.perStage.sh"
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
  
      submit "calibListCPass0" 1 1 "$LASTJOB" "${alirootEnv} ${self}"  PrintValues calibfile ${commonOutputPath}/meta/cpass0.calib.run${runNumber}.list "${commonOutputPath}/meta/cpass0.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="calibListCPass0"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit ${JOBID2} 1 1 "${LASTJOB}" "${alirootEnv} ${self}" MergeCPass0 ${targetDirectory} ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${commonOutputPath}/meta/cpass0.calib.run${runNumber}.list "${extraOpts[@]}"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID2}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      if [[ -n ${generateMC} ]]; then
        submit "mrl${JOBpostfix}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" PrintValues sim ${commonOutputPath}/meta/sim.run${runNumber}.list ${commonOutputPath}/meta/cpass0.job*.run${runNumber}.done
        # ---| set last job id used to submit dependency jobs |-------------------
        LASTJOB="mrl${JOBpostfix}"
        # treat the slurm case which uses the id number not the name
        # lastJobID is set in 'submit'
        if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
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
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID4wait}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
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
                        "${configPath}/mergeMakeOCDB.byComponent.perStage.sh"
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

      submit "calibListCPass1" 1 1 "$LASTJOB" "${alirootEnv} ${self}"  PrintValues calibfile ${commonOutputPath}/meta/cpass1.calib.run${runNumber}.list "${commonOutputPath}/meta/cpass1.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="calibListCPass1"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "qaListCPass1" 1 1 "$LASTJOB" "${alirootEnv} ${self}" PrintValues qafile ${commonOutputPath}/meta/cpass1.QA.run${runNumber}.lastMergingStage.txt.list "${commonOutputPath}/meta/cpass1.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="qaListCPass1"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "filteredListCPass1" 1 1 "$LASTJOB" "${alirootEnv} ${self}" PrintValues filteredTree ${commonOutputPath}/meta/cpass1.filtered.run${runNumber}.lastMergingStage.txt.list "${commonOutputPath}/meta/cpass1.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="filteredListCPass1"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "${JOBID5}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" MergeCPass1 ${targetDirectory} ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${commonOutputPath}/meta/cpass1.calib.run${runNumber}.list qaFilesToMerge=${commonOutputPath}/meta/cpass1.QA.run${runNumber}.lastMergingStage.txt.list filteredFilesToMerge=${commonOutputPath}/meta/cpass1.filtered.run${runNumber}.list "${extraOpts[@]}"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID5}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
      echo
    fi


    ################################################################################
    ################################################################################
    # run the CPass2 if requested

    if [ ${runCPass2reco} -eq 1 ]; then

      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass2"
      rm -f ${commonOutputPath}/meta/cpass2.job*.run${runNumber}.done

      # safety feature: if we are re-running for any reason we want to delete the previous output first.
      [[ -d ${targetDirectory} ]] && rm -rf ${targetDirectory}/* && echo "removed old output at ${targetDirectory}/*"

      echo
      echo "starting CPass2... for run ${runNumber}"
      echo

      # create directory and copy all files that are needed
      mkdir -p ${targetDirectory}
      
      filesCPass2=( 
                    "${configPath}/runPPass_pp.sh"
                    "${configPath}/runPPass_pbpb.sh"
                    "${configPath}/rec.C"
                    "${configPath}/runCalibTrain.C"
                    "${configPath}/QAtrain_duo.C"
                    "${configPath}/mergeQAgroups.C"
                    "${configPath}/localOCDBaccessConfig.C"
                    "${configPath}/OCDB.root"
      )
      for file in ${filesCPass2[*]}; do
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

      submit ${JOBID6} 1 ${nFiles} "${LASTJOB}" "${alirootEnv} ${self}" CPass2 ${targetDirectory} ${localInputList} ${nEvents} ${currentDefaultOCDB} ${configFile} ${runNumber} -1 "${extraOpts[@]}"

      ################################################################################
      ## submit a monitoring job that will run until a certain number of jobs are done with reconstruction
      submit "${JOBID6wait}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" WaitForOutput ${commonOutputPath} "meta/cpass2.job\*.run${runNumber}.done" ${nFilesToWaitFor} ${maxSecondsToWait}
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID6wait}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
      ################################################################################

      echo
    fi #end running CPass2

    ################################################################################
    # submit merging of CPass2, depends on the reconstruction
    if [ ${runCPass2MergeMakeOCDB} -eq 1 ]; then

      echo
      echo "submit CPass2 merging for run ${runNumber}"
      echo

      targetDirectory="${commonOutputPath}/${year}/${period}/000${runNumber}/cpass2"
      rm -f ${commonOutputPath}/meta/merge.cpass2.run${runNumber}.done
      mkdir -p ${targetDirectory}

      # copy files 
      filesMergeCPass2=(
                        "${configPath}/OCDB.root"
                        "${configPath}/localOCDBaccessConfig.C"
                        "${configPath}/mergeMakeOCDB.byComponent.perStage.sh"
                        "${configPath}/mergeMakeOCDB.byComponent.sh"
                        "${configPath}/mergeByComponent.C"
                        "${configPath}/makeOCDB.C"
                        "${configPath}/merge.C"
                        "${configPath}/mergeMakeOCDB.sh"
                        "${configPath}/QAtrain_duo.C"
			"${configPath}/mergeQAgroups.C"
      )
      for file in ${filesMergeCPass2[*]}; do
        [[ -f ${file} ]] && echo "copying ${file}" && cp -f ${file} ${commonOutputPath}
      done

      submit "calibListCPass2" 1 1 "$LASTJOB" "${alirootEnv} ${self}"  PrintValues calibfile ${commonOutputPath}/meta/cpass2.calib.run${runNumber}.list "${commonOutputPath}/meta/cpass2.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="calibListCPass2"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "qaListCPass2" 1 1 "$LASTJOB" "${alirootEnv} ${self}" PrintValues qafile ${commonOutputPath}/meta/cpass2.QA.run${runNumber}.lastMergingStage.txt.list "${commonOutputPath}/meta/cpass2.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="qaListCPass2"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "filteredListCPass2" 1 1 "$LASTJOB" "${alirootEnv} ${self}" PrintValues filteredTree ${commonOutputPath}/meta/cpass2.filtered.run${runNumber}.lastMergingStage.txt.list "${commonOutputPath}/meta/cpass2.job\*.run${runNumber}.done"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB="filteredListCPass2"
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

      submit "${JOBID7}" 1 1 "${LASTJOB}" "${alirootEnv} ${self}" MergeCPass2 ${targetDirectory} ${currentDefaultOCDB} ${configFile} ${runNumber} calibrationFilesToMerge=${commonOutputPath}/meta/cpass2.calib.run${runNumber}.list qaFilesToMerge=${commonOutputPath}/meta/cpass2.QA.run${runNumber}.lastMergingStage.txt.list filteredFilesToMerge=${commonOutputPath}/meta/cpass2.filtered.run${runNumber}.list "${extraOpts[@]}"
      # ---| set last job id used to submit dependency jobs |-------------------
      LASTJOB=${JOBID7}
      # treat the slurm case which uses the id number not the name
      # lastJobID is set in 'submit'
      if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
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
  ### final summary
  #################################################################################
  #if [ ${runESDfiltering} -eq 1 ]; then
  #  submit "${JOBID7wait}" 1 1 "${LASTJOB}" "${self}" WaitForOutput ${commonOutputPath} "meta/filtering.cpass2.run/\*.done" "${#listOfRuns[@]}" ${maxSecondsToWait}
  #else
  # TODO: most probably the file for the job 'meta/merge.cpass2.run\*.done' need to be adapted depending on the reconstrucion settings
    submit "${JOBID7wait}" 1 1 "${LASTJOB}" "${self}" WaitForOutput ${commonOutputPath} "meta/merge.cpass2.run\*.done" ${#listOfRuns[@]} ${maxSecondsToWait}
  #fi
  # ---| set last job id used to submit dependency jobs |-------------------
  LASTJOB=${JOBID7wait}
  # treat the slurm case which uses the id number not the name
  # lastJobID is set in 'submit'
  if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi

  #################################################################################
  echo
  echo "submit make a summary"
  echo

  [[ -z ${alirootEnvQA} ]] && alirootEnvQA=$(encSpaces "${alirootEnv}")
  submit "${JOBID8}" 1 1 "${LASTJOB}" "${alirootEnvQA} ${self}" MakeSummary ${configFile} "commonOutputPath=${commonOutputPath}"
  # ---| set last job id used to submit dependency jobs |-------------------
  LASTJOB=${JOBID8}
  # treat the slurm case which uses the id number not the name
  # lastJobID is set in 'submit'
  if [ -n "$lastJobID" ]; then LASTJOB=$lastJobID; fi
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

goSummaryToHTML() {(
set -e
# Header
cat <<EOF
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Release validation summary</title>
<style>
/* Body font */
body { font: 14px sans-serif; }

/* Links */
a { text-decoration: none;
    color: #CC444B; }

/* All tables */
table.report td { border-right: 2px solid #AAA; }
table.report tr > td:last-of-type { border-right: 0; }
table.report th { color: white; }
table.report th,
table.report td { padding: 7px 10px 7px 10px; }
table.report { border-collapse: collapse; }

/* All tables: alternate rows */
table.report tr:nth-of-type(odd) { background-color: white; }
table.report tr:nth-of-type(even) { background-color: #F0F2EF; }

/* Summary tables: align numbers */
table.summary td:nth-of-type(5) { text-align: right; }
table.fullsummary td:nth-of-type(3) { text-align: right; }

/* File in full summary */
tr.emph td { font-weight: bold;
             text-align: left !important; }

/* Table 1 */
div.group:nth-of-type(1) h2,
div.group:nth-of-type(1) h3 { color: #F63E02; }
div.group:nth-of-type(1) table.report th { background: #F63E02; }
div.group:nth-of-type(1) table.report td { border-right-color: #F63E02; }
div.group:nth-of-type(1) table.report tr.emph { color: #F63E02; }

/* Table 2 */
div.group:nth-of-type(2) h2,
div.group:nth-of-type(2) h3 { color: #FF6201; }
div.group:nth-of-type(2) table.report th { background: #FF6201; }
div.group:nth-of-type(2) table.report td { border-right-color: #FF6201; }
div.group:nth-of-type(2) table.report tr.emph { color: #FF6201; }

/* Table 3 */
div.group:nth-of-type(3) h2,
div.group:nth-of-type(3) h3 { color: #FAA300; }
div.group:nth-of-type(3) table.report th { background: #FAA300; }
div.group:nth-of-type(3) table.report td { border-right-color: #FAA300; }
div.group:nth-of-type(1) table.report tr.emph { color: #FAA300; }

/* Table 4 */
div.group:nth-of-type(4) h2,
div.group:nth-of-type(4) h3 { color: #F7B538; }
div.group:nth-of-type(4) table.report th { background: #F7B538; }
div.group:nth-of-type(4) table.report td { border-right-color: #F7B538; }
div.group:nth-of-type(1) table.report tr.emph { color: #F7B538; }

</style>
</head>
<body>
EOF

# Error summary
cat <<EOF
<div class="group">
<h2>Error summary</h2>
<table class="report summary">
<tr><th>Error type</th><th>Run</th><th>Pass</th><th>Files</th><th>Count</th></tr>
EOF
while read F; do
  F=$( echo $F | sed -e 's!\(.*\) *: run \([0-9]\+\) \([^ ]\+\) *\([^(]\+\) *( *\([0-9]\+\).*$!\1|\2|\3|\4|\5!g' )
  TYPE=$(echo $F|cut -d\| -f1)
  RUN=$(echo $F|cut -d\| -f2)
  PASS=$(echo $F|cut -d\| -f3|sed -e 's/cpass2/ppass/')
  FILES=$(echo $F|cut -d\| -f4|xargs echo)
  COUNT=$(echo $F|cut -d\| -f5)

  echo "<tr><td>$TYPE</td><td>$RUN</td><td>$PASS</td><td>$FILES</td><td>$COUNT</td></tr>"
done < <(awk '/error summary:/,/detailed summary:/' summary.log | grep -v '======' | grep ': run')
printf "</table>\n</div>\n"

# Detailed summary
echo "<div class=\"group\">"
echo "<h2>Detailed summary</h2>"
RUN=
PASS=
COUNT=
while read F; do
  [[ $F == run* ]] && { RUN=$(echo $F|cut -d' ' -f2); PASS=; continue; }
  if [[ $F == *cpass* ]]; then
    if [[ ! $PASS ]]; then
      [[ $COUNT ]] && echo "</table>"
      cat <<EOF
<h3>Run $RUN</h3>
<table class="report fullsummary">
<tr><th>Pass</th><th>File</th><th>Count</th></tr>
EOF
    fi
    PASS=$(echo $F|sed -e 's/[^a-z0-9 ]//g')
  elif [[ $F == '->ERRORS'* ]]; then
    F=$(echo $F|cut -d: -f2-)
    FILE=$(echo $F|cut -d' ' -f1)
    OK=$(echo $F|cut -d' ' -f2|sed -e 's/OK://g')
    BAD=$(echo $F|cut -d' ' -f3|sed -e 's/BAD://g')
    [[ $OK ]] || OK=0
    [[ $BAD ]] || BAD=0
    echo "<tr class=\"emph\"><td>$FILE</td><td>$OK OK</td><td>$BAD bad</td></tr>"
  else
    COUNT=$(echo $F|cut -d' ' -f1)
    FILE=$(echo $F|cut -d' ' -f3-)
    echo "<tr><td>$PASS</td><td>$FILE</td><td>$COUNT</td></tr>"
  fi
done < <(awk '/detailed summary:/,/list of bad logs:/' summary.log | grep -v '======' | grep -v -- '___' | sed -e '/^$/d')
[[ $COUNT ]] && echo "</table>"
echo "</div>"

# Bad logs
cat <<EOF
<div class="group">
<h2>Logs with errors</h2>
<table class="report badfiles">
<tr><th>Log file</th><th>Error type</th></tr>
EOF
while read F; do
  FILE=$(echo $F|cut -d' ' -f1|sed -e 's!^.*/\(AliPhysics-[^/]\+\)/!!g')
  ERROR=$(echo $F|cut -d' ' -f2-)
  FILE="<a href=\"$FILE\">$FILE</a>"
  echo "<tr><td>$FILE</td><td>$ERROR</td></tr>"
done < <(awk '/bad logs:/,/core files:/' summary.log | grep -v '======' | grep ://)
printf "</table>\n</div>\n"

# Core files
cat <<EOF
<div class="group">
<h2>Core files</h2>
<table class="report badfiles">
<tr><th>Core file</th></tr>
EOF
while read F; do
  FILE=$(echo $F|cut -d' ' -f1|sed -e 's!^.*/\(AliPhysics-[^/]\+\)/!!g')
  FILE="<a href=\"$FILE\">$FILE</a>"
  echo "<tr><td>$FILE</td></tr>"
done < <(awk '/core files:/,/^$/' summary.log | grep -v '======' | grep ://)
printf "</table>\n</div>\n"

# Footer
cat <<EOF
</body>
</html>
EOF
set +e
)}

goMakeSummary()
(
  # All the final stuff goes in here for ease of use:
  #  - summary logs
  #  - qa plot making
  #  - final file lists
  # Runs in current dir - in makeflow mode it can run LOCAL, then the QA plots and summaries
  # will appear in the submission dir.

  #some defaults:
  log=summary_full.log
  productionID=qa
  alilog_info "[BEGIN] goMakeSummary() with following parameters $*"
  configFile=$1
  shift 1
  extraOpts=("$@")

  parseConfig configFile=$configFile "${extraOpts[@]}" || return 1
  batchWorkingDirectory=$PWD

  logTmp="${batchWorkingDirectory}/${log}"

  [[ -f "$alirootSource" && -z "$ALICE_ROOT" ]] && source "$alirootSource"
  [[ ! -f "$configFile" ]] && echo "no config file ${configFile}!" && return 1

  # If we use workqueue all the files are in PWD, on batch systems they are in the meta subdir.
  [[ -z ${commonOutputPath} ]] && commonOutputPath=${batchWorkingDirectory}
  metadir=${commonOutputPath}/meta
  [[ ! -d $metadir ]] && metadir="$PWD"

  # Redirect all output to both file and console. Save fd 3 to restore later.
  exec 3>&1
  exec &> >(tee ${logTmp})

  # Take a snapshot of the current directory. Files already here at this point needn't be copied to
  # the destination, as they were transferred here by Makeflow/Work Queue.
  dirSnapshotExclusion=$(for f in *; do echo -n "$f|"; done)
  dirSnapshotExclusion="${log}"  # Append logfile. We'll copy it separately. Can't copy if in use.

  echo "[goMakeSummary] Exclusion list: those files are already here and will not be copied to destination: $dirSnapshotExclusion"
  listDir "$PWD"

  # Summarize the global stuff
  echo "env script: ${alirootSource} ${alirootEnv}"
  echo "ALICE_ROOT=${ALICE_ROOT}"
  echo "commonOutputPath=${commonOutputPath}"

  # Report simple summary (as opposed to a bloated log) to a separated file.
  #  * Full log will be on summary_full.log.
  #  * Simplified log will be on summary.log.
  rm -f summary.log
  goSummarizeMetaFiles "$metadir" | tee -a summary.log
  goSummaryToHTML > summary.html

  if [[ $simplifiedSummary == 1 ]]; then
    maxCopyTries=1 remoteCpTimeout=120 xCopy -f -d $PWD $commonOutputPath/running_time
    { [[ -e running_time ]] && cat running_time || echo "No running time info available."; } | tee -a summary.log
    alilog_info "[END] goMakeSummary() (stopping after simplified summary) with following parameters $*"
    exec 1>&3 3>&-
    exec 2>&1
    xCopy -f -d $commonOutputPath summary.log summary_full.log summary.html
    goCheckSummary
    return $?
  fi

  #make file lists - QA, trending, stacktraces, filtering and calibration
  ### wait for the merging of all runs to be over ###
  goPrintValues qafile remote.cpass1.qa.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues calibfile remote.calib.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues trendingfile remote.trending.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues filteredTree remote.cpass1.filtering.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues dcsTree remote.cpass0.dcsTree.list "$metadir"/merge.cpass0.run*.done &>/dev/null
  goPrintValues dcsTree remote.cpass1.dcsTree.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues stacktrace remote.cpass0.stacktrace.list "$metadir"/*cpass0*done &>/dev/null
  goPrintValues stacktrace remote.cpass1.stacktrace.list "$metadir"/*cpass1*done &>/dev/null
  goPrintValues stacktrace remote.cpass2.stacktrace.list "$metadir"/*cpass2*done &>/dev/null
  goPrintValues esd remote.cpass1.esd.list "$metadir"/*cpass1*done &> /dev/null
  goPrintValues syswatchRec remote.cpass0.syswatch.rec.list "$metadir"/merge.cpass0.run*.done &>/dev/null
  goPrintValues syswatchCalib remote.cpass0.syswatch.calib.list "$metadir"/merge.cpass0.run*.done &>/dev/null
  goPrintValues syswatchRec remote.cpass1.syswatch.rec.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues syswatchCalib remote.cpass1.syswatch.calib.list "$metadir"/merge.cpass1.run*.done &>/dev/null
  goPrintValues syswatchRec remote.cpass2.syswatch.rec.list "$metadir"/merge.cpass2.run*.done &>/dev/null
  goPrintValues qafile remote.cpass2.qa.list "$metadir"/merge.cpass2.run*.done &>/dev/null

  listDir "$PWD" "after goPrintValues"

  # Copy all the files to a local dir tree.
  for remoteList in remote.cpass1.qa.list \
                    remote.calib.list \
                    remote.trending.list \
                    remote.cpass1.filtering.list \
                    remote.cpass0.dcsTree.list \
                    remote.cpass1.dcsTree.list \
                    remote.cpass0.stacktrace.list \
                    remote.cpass1.stacktrace.list \
                    remote.cpass2.stacktrace.list \
                    remote.cpass2.qa.list \
                    remote.cpass{0,1,2}.syswatch.rec.list; do
    localList=${remoteList#remote.}
    rm -f "$localList" && touch "$localList"
    while read sourceFile; do
      destinationFile="${PWD}/${sourceFile#${commonOutputPath}}"
      copyFileFromRemote "$sourceFile" "$(dirname "${destinationFile}")" && \
        echo "$destinationFile" >> "$localList"
    done < <(cat "$remoteList")
  done

  listDir "$PWD" "after copying remote lists locally"

  #summarize the stacktraces
  stackTraceTree @stacktraces.cpass0.list > stacktrace_cpass0.tree
  stackTraceTree @stacktraces.cpass1.list > stacktrace_cpass1.tree
  stackTraceTree @stacktraces.cpass2.list > stacktrace_cpass2.tree

  #merge syslogs
  mergeSysLogs syslog.rec.cpass0.tree @cpass0.syswatch.rec.list
  mergeSysLogs syslog.calib.cpass0.tree @cpass0.syswatch.calib.list
  mergeSysLogs syslog.rec.cpass1.tree @cpass1.syswatch.rec.list
  mergeSysLogs syslog.calib.cpass1.tree @cpass1.syswatch.calib.list
  mergeSysLogs syslog.rec.cpass2.tree @cpass2.syswatch.rec.list

  #merge trending
  rm -f trending.root
  goMerge trending.list trending.root ${configFile} "${extraOpts[@]}" &> mergeTrending.log

  listDir "$PWD" "after merging trending"

  printExec goMakeSummaryTree "${metadir}" 0
  printExec goMakeSummaryTree "${metadir}" 1

  listDir "$PWD" "after making summary tree"

  if [[ -z "$pretend" ]]; then
    goCreateQAplots "${PWD}/cpass1.qa.list" "${productionID}" "QAplots_CPass1" "${configFile}" "${extraOpts[@]}" filteringList="${PWD}/cpass1.filtering.list" &>createQAplots.cpass1.log
    goCreateQAplots "${PWD}/cpass2.qa.list" "${productionID}" "QAplots_CPass2" "${configFile}" "${extraOpts[@]}" filteringList="${PWD}/cpass2.filtering.list" &>createQAplots.cpass2.log
  else
    mkdir QAplots_CPass2
    touch QAplots_CPass2/log
    mkdir QAplots_CPass1
    touch QAplots_CPass1/log
  fi

  listDir "$PWD" "after creation of QA plots"

  #make a merged summary tree out of the QA trending, dcs trees and log summary trees
  goMakeMergedSummaryTree

  listDir "$PWD" "after making merged summary tree"

  #if set, email the summary
  [[ -n ${MAILTO} ]] && cat ${logTmp} | mail -s "benchmark ${productionID} done" ${MAILTO}

  # Copy all, recursively, with the exception of snapshotted files already present when this
  # function was called.
  [[ "$dirSnapshotExclusion" == '' ]] \
    && xCopy -d $commonOutputPath/ . \
    || xCopy -d $commonOutputPath/ ./!($dirSnapshotExclusion)

  # Copy stdout to destination.
  echo "Copying ${log}. NOTE: this is the last bit you will see in the log!"
  exec 1>&3 3>&-
  exec 2>&1
  xCopy -f -d $commonOutputPath/ summary.log summary_full.log

  alilog_info "[END] goMakeSummary() with following parameters $*"
  return 0
)

goSummarizeMetaFiles()
{
  #summarize the meta files in current dir (unless specified differently)
  find ${1-"."} -name "*.done" -exec echo donefile {} \; -exec cat {} \; | \
  gawk '
    BEGIN {
    }

    /^donefile /  {
      donefile=$2
      sub(/^.*\//,"",donefile)
      match($2,/.pass[0-9]?/,tmparr)
      pass=tmparr[0]
      match($2,/run([0-9]*)/,tmparr)
      runNumber=tmparr[1]
      runs[runNumber]=1;
      passes[pass]=1;
    }

    #logfiles
    /OK/ || /BAD/ {
      logFile=$1
      sub(/^.*\//,"",logFile)

      if (donefile ~ /merge\./) {
        mergeLogsFiles[logFile]=1
        mergeLogs[runNumber,pass,logFile]++
        if ($0 ~ /OK/) mergeLogsOK[runNumber,pass,logFile]++
        else if ($0 ~ /BAD/) {
          mergeLogsBAD[runNumber,pass,logFile]++
          listOfBadLogs[nBadLogs++]=$0
        }
      }

      else if (donefile ~ /\.job/) {
        jobLogsFiles[logFile]=1
        jobLogs[runNumber,pass,logFile]++
        if ($0 ~ /OK/) jobLogsOK[runNumber,pass,logFile]++
        else if ($0 ~ /BAD/) {
          jobLogsBAD[runNumber,pass,logFile]++
          listOfBadLogs[nBadLogs++]=$0
        }
      }

      next
    }

    #output files
    {
      fileType=$1
      file=$2
      if (donefile ~ /merge\./) {
        mergeFileTypes[fileType]=1
        outputFilesMerge[runNumber,pass,fileType]++
      }
      else if (donefile ~ /\.job/) {
        jobFileTypes[fileType]=1
        outputFilesJobs[runNumber,pass,fileType]++
      }
      if (fileType == "core") coreFiles[nCoreFiles++]=file
    }

    END {
      print "===== error summary: ================================================================="
      for (run in runs ) {
        for (pass in passes ) {

          for (logFile in jobLogsFiles) {
            if (jobLogsBAD[run,pass,logFile]>0) {
              print "ERROR      : run "run" "pass" "logFile" ( "jobLogsBAD[run,pass,logFile]" failures )"
            }
          }

          for (logFile in mergeLogsFiles) {
            if (mergeLogsBAD[run,pass,logFile]>0) {
              print "ERROR merge: run "run" "pass" "logFile" ( "mergeLogsBAD[run,pass,logFile]" failures )"
            }
          }

          if (outputFilesJobs[run,pass,"core"]>0) {
            print "CORE       : run "run" "pass" ( "outputFilesJobs[run,pass,"core"]" core files! )"
          }
          if (outputFilesMerge[run,pass,"core"]>0) {
            print "CORE       : run "run" merge "pass" ( "outputFilesMerge[run,pass,"core"]" core files! )"
          }

        }
      }

      print ""
      print "===== detailed summary: ================================================================="
      for (run in runs ) {
        print"________________________________________"
        print "run "run

        for (pass in passes ) {
          print "=="pass"=="

          filesToPrint["esd"]=0
          filesToPrint["qafile"]=0
          filesToPrint["ocdbTarball"]=0
          filesToPrint["aod"]=0
          filesToPrint["core"]=0
          for (outputFile in filesToPrint) {
            if (outputFilesJobs[run,pass,outputFile])  print "      "outputFilesJobs[run,pass,outputFile]" X "outputFile
          }
          print "  "pass" merge:"
          for (outputFile in filesToPrint) {
            if (outputFilesMerge[run,pass,outputFile])  print "      "outputFilesMerge[run,pass,outputFile]" X "outputFile
          }
          print ""
          for (logFile in jobLogsFiles) {
            if (jobLogsBAD[run,pass,logFile]>0) {
              print "    ->ERRORS: "logFile" OK:"jobLogsOK[run,pass,logFile]" BAD:"jobLogsBAD[run,pass,logFile]
            }
          }

          for (logFile in mergeLogsFiles) {
            if (mergeLogsBAD[run,pass,logFile]>0) {
              print "    ->ERRORS merge: "logFile" OK:"mergeLogsOK[run,pass,logFile]" BAD:"mergeLogsBAD[run,pass,logFile]
            }
          }

          print ""
        }
      }

      if (nBadLogs)
      {
        print ""
        print "===== list of bad logs: ================================================================="
        for (i in listOfBadLogs) {
          print listOfBadLogs[i]
        }
      }
      if (nCoreFiles)
      {
        print ""
        print "===== list of core files: ================================================================="
        for (i in coreFiles) {
          print coreFiles[i]
        }
      }
    }'
}

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
