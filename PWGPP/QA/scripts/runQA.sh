#!/usr/bin/env bash
# process QA output into plots and trending
# run without arguments for examples
# origin: Mikolaj Krzewicki, mkrzewic@cern.ch
#
if [ ${BASH_VERSINFO} -lt 4 ]; then
  echo "bash version >= 4 needed, you have ${BASH_VERSION}, exiting..."
  exit 1
fi

#run in proper mode depending on the selection
for scr in "utilities.sh" "alilog4bash.sh"; do
  if [[ -e $scr ]]; then
    echo "Sourcing $scr from current directory"
    source $scr false
  elif [[ -e $ALICE_ROOT/libexec/$scr ]]; then
    echo "Sourcing $scr from AliRoot"
    source $ALICE_ROOT/libexec/$scr false
  else
    exit 1
  fi
done

[[ -z $ALICE_ROOT ]] && echo "ALICE_ROOT not defined" && exit 1

##############################
#default values:
##############################
#config file
configFile=""
#where to search for qa files
inputList=file.list
#working directory
workingDirectory="${PWD}"
#where to place the final qa plots
#outputDirectory="/afs/cern.ch/work/a/aliqa%det/www/"
outputDirectory="${workingDirectory}/%DET"
#filter out detector option
excludeDetectors="EXAMPLE"
#logs
logDirectory=${workingDirectory}/logs
#OCDB storage
ocdbStorage="raw://"
#arbitrary id
productionID="default"
#send email as:
MAILFROM=""
#email to
MAILTO=""
MAILTO_TPC=""
#mail a short status message to "MAILTO" when ready
MAILshortSummary=""
#attach debug info
MAILdebugInfo=1
MAILfullProductionLog=1
MAILcompressLogs=""

main()
{
  if [[ -z $1 ]]; then
    echo "Usage: "
    echo "  ${0##*/} option=value [option=value]"
    echo "  at least inputList should be specified, or configFile containing it:"
    echo "  ${0##*/} inputList=file.list"
    echo "  options override config file (if any), e.g.:"
    echo "  ${0##*/} configFile=runQA.config inputList=file.list outputDirectory=%det"
    echo "some expert options"
    echo "  inputListHighPtTrees=file.list - external list of filtered trees, requires inputList to be set"
    echo "  includeDetectors=TPC,V0,MU - only process those"
    echo "  excludeDetectors=EVS,TPC - skip processing of those"
    echo "  - see example config file for more"
    return 1
  fi

  if ! parseConfig "$@"; then
    return 1
  fi

  ocdbregex='raw://'
  if [[ ${ocdbStorage} =~ ${ocdbregex} ]]; then
    alien-token-init ${alienUserName}
    #this is a hack! alien-token init seems not enough
    #but the gclient_env script messes up the LD_LIBRARY_PATH
    while read x; do
      eval ${x};
    done < <(grep -v "LD_LIBRARY_PATH" /tmp/gclient_env_${UID})
  fi

  updateQA "$@"
  return 0
}

updateQA()
{
  umask 0002
  parseConfig "$@"

  #be paranoid and make some full paths
  [[ ! -f ${inputList} ]] && echo "no input list: ${inputList}" && return 1
  inputList=$(get_realpath ${inputList})
  mkdir -p ${workingDirectory}
  #this is a trick to get the full path of workingDirectory
  #(on a mac 'readlink -f' does not work...)
  workingDirectory=$(workingDirectory=${workingDirectory%/}; cd ${workingDirectory%/*}; echo "${PWD}/${workingDirectory##*/}")
  if [[ ! -d ${workingDirectory} ]]; then
    echo "working dir $workingDirectory does not exist and cannot be created"
    return 1
  fi
  cd ${workingDirectory}

  echo JOB config:
  echo inputList=$inputList
  echo outputDirectory=$outputDirectory
  echo

  dateString=$(date +%Y-%m-%d-%H-%M-%S-%N)
  echo "Start time QA process: $dateString"
  uniquePID=$(createUniquePID "${productionID}")

  #logging
  mkdir -p $logDirectory
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && return 1
  logFile="$logDirectory/${0##*/}.${dateString}.log"
  logFileShort="$logDirectory/${0##*/}.${dateString}.short.log"
  touch ${logFile}
  [[ ! -f ${logFile} ]] && echo "cannot write logFile $logFile" && return 1
  echo "logFile = $logFile"

  #check lock
  lockFile=${workingDirectory}/runQA.lock
  [[ -f ${lockFile} ]] && echo "lock ${lockFile} exists!" | tee ${logFile} && return 1
  touch ${lockFile}
  [[ ! -f ${lockFile} ]] && echo "cannot lock $lockFile" | tee ${logFile} && return 1

  exec &>${logFile}

  ################################################################
  #ze detector loop
  declare -A arrLogSummary
  declare -A arrDetectorStatus
  for detectorScript in $ALICE_PHYSICS/PWGPP/QA/detectorQAscripts/*; do
    echo
    echo "##############################################"
    echo $(date)
    detectorStatus="OK"
    [[ ! ${detectorScript} =~ .*\.sh$ ]] && continue
    detector=${detectorScript%.sh}
    detector=${detector##*/}
    #by default we expect the container in the QA root file to de named like
    #the detector
    detectorQAcontainerName=${detector}

    #reset to the default value the options which can be overridden
    #by the detectors in their scripts
    MAILdebugInfo=${parseConfig__ORIGINAL__MAILdebugInfo}
    MAILcompressLogs=${parseConfig__ORIGINAL__MAILcompressLogs}
    MAILfullProductionLog=${parseConfig__ORIGINAL__MAILfullProductionLog}

    #skip if excluded
    if [[ "${excludeDetectors}" =~ ${detector} ]]; then
      echo "${detector} is excluded in config, skipping..."
      continue
    fi

    #if includeDetectors set, only process thoe detectors specified there
    if [[ -n ${includeDetectors} && ! "${includeDetectors}" =~ ${detector} ]]; then
      echo "${detector} not included in includeDetectors, skipping..."
      continue
    fi

    logSummary=${logDirectory}/summary-${detector}-${dateString}.log
    hostInfo >> ${logSummary}
    outputDir=$(substituteDetectorName ${detector} ${outputDirectory})
    tmpDetectorRunDir=${workingDirectory}/tmpQAtmpRunDir${detector}-${dateString}
    if ! mkdir -p ${tmpDetectorRunDir}; then
      echo "cannot create the temp dir $tmpDetectorRunDir"
      continue
    fi
    #at this point we can store the detector name and the logSummary location
    #in arrLogSummary
    arrLogSummary[${detector}]=${logSummary}

    cd ${tmpDetectorRunDir}

    tmpPrefix=${tmpDetectorRunDir}/${outputDir}
    echo "running QA for ${detector}"
    echo "  outputDir=$outputDir"
    echo "  tmpPrefix=$tmpPrefix"

    #source the detector script
    #unset the detector functions from previous iterations (detectors)
    unset -f runLevelQA
    unset -f runLevelAodQA
    unset -f runLevelQAouter
    unset -f periodLevelQA
    unset -f runLevelEventStatQA
    unset -f runLevelHighPtTreeQA
    unset -f periodLevelHighPtTreeQA
    source ${detectorScript}

    #################################################################
    #produce the QA and trending tree for each file (run)
    unset arrOfTouchedProductions
    declare -A arrOfTouchedProductions
    while read inputFile; do
      echo
      echo $(date)
      echo "INFO Input file is:" $inputFile
      #first check if input file exists
      [[ ! -f ${inputFile%\#*} ]] && echo "file ${inputFile%\#*} not accessible" && continue

      if ! guessRunData ${inputFile}; then
        echo "could not guess run data from ${inputFile}"
        continue
      fi
      echo "anchorYear for ${originalPeriod} is: ${anchorYear}"

      if [[ ${dataType} =~ "sim" ]]; then
        tmpProductionDir=${tmpPrefix}/${dataType}/${year}/${originalPeriod}/${pass}
      else
        tmpProductionDir=${tmpPrefix}/${dataType}/${year}/${period}/${pass}
      fi

      tmpRunDir=${tmpProductionDir}/000${runNumber}
      mkdir -p ${tmpRunDir}

      cd ${tmpRunDir}

      #check what kind of input file we have, default is a zip archive
      #set the inputs accordingly
      qaFile=""
      qaFileAod=""
      qaFileOuter=""
      highPtTree=""
      eventStatFile=""
      eventStatFileOuter=""
      #it is possible we get the highPt trees from somewhere else
      #search the list of high pt trees for the proper run number
      if [[ -n ${inputListHighPtTrees} ]]; then
        highPtTree=$(egrep -m1 ${runNumber} ${inputListHighPtTrees})
        echo "loaded the highPtTree ${highPtTree} from external file ${inputListHighPtTrees}"
      fi
      #if we are explicit about the input file this takes precedence
      #over earlier additions
      [[ "${inputFile}" =~ QAresults.root$ ]] && qaFile=${inputFile}
      [[ "${inputFile}" =~ QAresults_AOD.root$ ]] && qaFileAod=${inputFile}
      [[ "${inputFile}" =~ QAresults_merged.root$ ]] && qaFile=${inputFile}
      [[ "${inputFile}" =~ QAresults_barrel.root$ ]] && qaFile=${inputFile}
      [[ "${inputFile}" =~ QAresults_outer.root$ ]] && qaFileOuter=${inputFile}
      [[ "${inputFile}" =~ FilterEvents_Trees.root$ ]] && highPtTree=${inputFile}
      [[ "${inputFile}" =~ event_stat.root$ ]] && eventStatFile=${inputFile}
      [[ "${inputFile}" =~ event_stat_barrel.root$ ]] && eventStatFile=${inputFile}
      [[ "${inputFile}" =~ event_stat_outer.root$ ]] && eventStatFileOuter=${inputFile}
      if [[ "${inputFile}" =~ \.zip$ ]]; then
        [[ -z ${qaFile} ]] && qaFile=${inputFile}
        [[ -z ${qaFileAod} ]] && qaFileAod=${inputFile}
        [[ -z ${qaFileOuter} ]] && qaFileOuter=${inputFile}
        [[ -z ${highPtTree} ]] && highPtTree=${inputFile}
        [[ -z ${eventStatFile} ]] && eventStatFile=${inputFile}
        [[ -z ${eventStatFileOuter} ]] && eventStatFileOuter=${inputFile}
      fi

      #if we have zip archives in the input, extract the proper file name
      #from the archive and append in a root-like fashion
      if [[ "$qaFile" =~ .*.zip$ ]]; then
        if unzip -l ${qaFile} | egrep "QAresults.root" &>/dev/null; then
          qaFile+="#QAresults.root"
        elif unzip -l ${qaFile} | egrep "QAresults_barrel.root" &>/dev/null; then
          qaFile+="#QAresults_barrel.root"
        else
          qaFile=""
        fi
      fi
      if [[ "$qaFileOuter" =~ .*.zip$ ]]; then
        if unzip -l ${qaFileOuter} | egrep "QAresults_outer.root" &>/dev/null; then
          qaFileOuter+="#QAresults_outer.root"
        else
          qaFileOuter=""
        fi
      fi
      if [[ "$highPtTree" =~ .*.zip$ ]]; then
        if unzip -l ${highPtTree} | egrep "FilterEvents_Trees.root" &>/dev/null; then
          highPtTree+="#FilterEvents_Trees.root"
        else
          highPtTree=""
        fi
      fi
      if [[ "${eventStatFile}" =~ .*.zip$ ]]; then
        if unzip -l ${eventStatFile} | egrep "event_stat.root" &>/dev/null; then
          eventStatFile+="#event_stat.root"
        elif unzip -l ${eventStatFile} | egrep "event_stat_barrel.root" &>/dev/null; then
          eventStatFile+="#event_stat_barrel.root"
        else
          eventStatFile=""
        fi
      fi
      if [[ "${eventStatFileOuter}" =~ .*.zip$ ]]; then
        if unzip -l ${eventStatFileOuter} | egrep "event_stat_outer.root" &>/dev/null; then
          eventStatFileOuter+="#event_stat.root"
        else
          eventStatFileOuter=""
        fi
      fi
     if [[ "${qaFileAod}" =~ .*.zip$ ]]; then
        if unzip -l ${qaFileAod} | egrep "QAresults_AOD.root" &>/dev/null; then
          qaFileAod+="#QAresults_AOD.root"
        else
          qaFileAod=""
        fi
      fi
  
      echo qaFile=$qaFile
      echo qaFileAod=$qaFileAod
      echo qaFileOuter=$qaFileOuter
      echo highPtTree=$highPtTree
      echo eventStatFile=$eventStatFile
      echo eventStatFileOuter=$eventStatFileOuter
      echo ocdbStorage=${ocdbStorage}
      echo

      #standard QA based on QAresults.root file (and variants)
      if [[ -n ${qaFile} && $(type -t runLevelQA) =~ "function" ]]; then
        echo running ${detector} runLevelQA for run ${runNumber} from ${qaFile}
        ( runLevelQA "${qaFile}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi

     #standard QA based on QAresults_AOD.root file (and variants)
      if [[ -n ${qaFileAod} && $(type -t runLevelAodQA) =~ "function" ]]; then
        echo running ${detector} runLevelAodQA for run ${runNumber} from ${qaFileAod}
        ( runLevelAodQA "${qaFileAod}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi



      #standard QA based on QAresults_outer.root file (there in cpass, with different triggers)
      if [[ -n ${qaFileOuter} && $(type -t runLevelQAouter) =~ "function" ]]; then
        echo running ${detector} runLevelQAouter for run ${runNumber} from ${qaFileOuter}
        ( runLevelQAouter "${qaFileOuter}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi
      #expert QA based on high pt trees
      if [[ -n ${highPtTree} && $(type -t runLevelHighPtTreeQA) =~ "function" ]]; then
        echo running ${detector} runLevelHighPtTreeQA for run ${runNumber} from ${highPtTree}
        ( runLevelHighPtTreeQA "${highPtTree}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi
      #event stat QA based on event_stat.root file
      if [[ -n ${eventStatFile} && $(type -t runLevelEventStatQA) =~ "function" ]]; then
        echo running ${detector} runLevelEventStatQA for run ${runNumber} from ${eventStatFile}
        ( runLevelEventStatQA "${eventStatFile}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi
      #event stat QA based on event_stat_outer.root file
      if [[ -n ${eventStatFileOuter} && $(type -t runLevelEventStatQAouter) =~ "function" ]]; then
        echo running ${detector} runLevelEventStatQAouter for run ${runNumber} from ${eventStatFileOuter}
        ( runLevelEventStatQAouter "${eventStatFileOuter}" ) &>> runLevelQA.log
        #cache the touched production + an example file to guarantee consistent run data parsing
        arrOfTouchedProductions[${tmpProductionDir}]="${inputFile%\#*}"
      fi

      #perform some default actions:
      #if trending.root not created, create a default one, if anything was actually processed
      if [[ ! -f trending.root && -n ${qaFile} ]]; then
        echo "trending.root not provided, falling back to:"
        echo "aliroot -b -q -l $ALICE_PHYSICS/PWGPP/macros/simpleTrending.C(\"${qaFile}\",${runNumber},\"${detectorQAcontainerName}\",\"trending.root\",\"trending\",\"recreate\") 2>&1 | tee -a runLevelQA.log > simpleTrending.log"
        aliroot -b -q -l "$ALICE_PHYSICS/PWGPP/macros/simpleTrending.C(\"${qaFile}\",${runNumber},\"${detectorQAcontainerName}\",\"trending.root\",\"trending\",\"recreate\")" 2>&1 | tee -a runLevelQA.log > simpleTrending.log
      fi
      if [[ ! -f trending.root ]]; then
        echo "trending.root not created"
      fi

      cd ${tmpDetectorRunDir}

    done < ${inputList}

    #################################################################
    #cache which productions were (re)done
    echo "list of processed productions:"
    echo "    ${!arrOfTouchedProductions[@]}"
    echo

    #################################################################
    #(re)do the merging/trending
    for tmpProductionDir in ${!arrOfTouchedProductions[@]}; do
      cd ${tmpProductionDir}
      echo
      echo "running period level stuff in ${tmpProductionDir}"
      echo $(date)

      productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}
      echo productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}

      mkdir -p ${productionDir}
      if [[ ! -d ${productionDir} ]]; then
        echo "cannot make productionDir $productionDir" && continue
      fi

      #move runs to final destination
      for dir in ${tmpProductionDir}/000*; do
        echo
        oldRunDir=${outputDir}/${dir#${tmpPrefix}}
        if ! guessRunData "${arrOfTouchedProductions[${tmpProductionDir}]}"; then
          echo "could not guess run data from ${arrOfTouchedProductions[${tmpProductionDir}]}"
          continue
        fi

        #before moving - VALIDATE!!!
        if ! validate ${dir}; then
          continue
        fi

        #moving a dir is an atomic operation, no locking necessary
        if [[ -d ${oldRunDir} ]]; then
          echo "removing old ${oldRunDir}"
          rm -rf ${oldRunDir}
        fi
        echo "moving new ${runNumber} to ${productionDir}"
        mv -f ${dir} ${productionDir}
      done

      #go to a temp dir to do the period level stuff in a completely clean dir
      tmpPeriodLevelQAdir="${tmpProductionDir}/periodLevelQA"
      echo
      echo tmpPeriodLevelQAdir="${tmpPeriodLevelQAdir}"
      if ! mkdir -p ${tmpPeriodLevelQAdir}; then continue; fi
      cd ${tmpPeriodLevelQAdir}

      #link the final list of per-run dirs here, just the dirs
      #to have a clean working directory
      unset linkedStuff
      declare -a linkedStuff
      for x in ${productionDir}/000*; do [[ -d $x ]] && ln -s $x && linkedStuff+=(${x##*/}); done

      #merge trending files if any
      if /bin/ls 000*/trending.root &>/dev/null; then # skip corrupted and empty trending files
        hadd -k trending.root 000*/trending.root &> periodLevelQA.log
      fi

      #run the period level trending/QA
      if [[ -f "trending.root" && $(type -t periodLevelQA) =~ "function" ]]; then
        echo running ${detector} periodLevelQA for production ${period}/${pass}
        ( periodLevelQA trending.root ) &>> periodLevelQA.log
      else
        echo "WARNING: not running ${detector} periodLevelQA for production ${period}/${pass}, no trending.root"
      fi

      if ! validate ${PWD}; then continue; fi

      #here we are validated so move the produced QA to the final place
      #clean up linked stuff first
      [[ -n ${linkedStuff[@]} ]] && rm ${linkedStuff[@]}
      periodLevelLock=${productionDir}/runQA.lock
      if [[ ! -f ${periodLevelLock} ]]; then
        #some of the output could be a directory, so handle that
        #TODO: maybe use rsync?
        #lock to avoid conflicts:
        echo "${uniquePID}" > ${periodLevelLock}
        for x in ${tmpPeriodLevelQAdir}/*; do
          if [[ -d ${x} ]]; then
            echo "removing ${productionDir}/${x##*/}"
            rm -rf ${productionDir}/${x##*/}
            echo "moving ${x} to ${productionDir}"
            mv ${x} ${productionDir}
          fi
          if [[ -f ${x} ]]; then
            echo "moving ${x} to ${productionDir}"
            mv -f ${x} ${productionDir}
          fi
        done
        rm -f ${periodLevelLock}
        #remove the temp dir
        rm -rf ${tmpPeriodLevelQAdir}
      else
        echo "ERROR: cannot move to destination"                     >> ${logSummary}
        echo "production dir ${productionDir} locked!"               >> ${logSummary}
        echo "check and maybe manually do:"                          >> ${logSummary}
        echo " rm ${periodLevelLock}"                                >> ${logSummary}
        echo " rsync -av ${tmpPeriodLevelQAdir}/ ${productionDir}/"  >> ${logSummary}
        detectorStatus="planB"
      fi

    done #end of merging/trending loop

    cd ${workingDirectory}

    #set the detectorStatus flag for the detector and proceed with cleanup
    arrDetectorStatus[${detector}]="${detectorStatus}"
    if [[ "${detectorStatus}" == "OK" ]]; then
      echo
      echo "detectorStatus=${detectorStatus}"
      echo removing ${tmpDetectorRunDir}
      rm -rf ${tmpDetectorRunDir}
    fi
  done #end of detector loop

  #make a run/log summary
  #make one big stacktrace tree for core files
  stackTraceArr=($(awk '/stacktrace\.log/ {print $1}' ${logDirectory}/summary-*-${dateString}.log))
  if [[ ${#stackTraceArr[@]} > 0 ]]; then
    stackTraceTree "${stackTraceArr[@]}" > ${logDirectory}/stacktrace-core-${dateString}.tree
  fi
  #sometimes core files are not available, use the root logs, which sometimes do have a stacktrace
  stackTraceArr=($(awk '/log[[:space:]]*BAD / {print $1}' ${logDirectory}/summary-*-${dateString}.log))
  if [[ ${#stackTraceArr[@]} > 0 ]]; then
    stackTraceTree "${stackTraceArr[@]}" > ${logDirectory}/stacktrace-log-${dateString}.tree
  fi
  #make stacktrace plots
  plotStackTraceTree ${logDirectory}/stacktrace-core-${dateString}.tree ${logDirectory}/stacktrace-core-${dateString}.png
  plotStackTraceTree ${logDirectory}/stacktrace-log-${dateString}.tree ${logDirectory}/stacktrace-log-${dateString}.png

  #process alarms, send emails etc, for each detector separately
  echo ""
  echo "log statistics:"
  for detector in "${!arrLogSummary[@]}"; do
    echo "${detector}"
    logSummary="${arrLogSummary[${detector}]}"
    printLogStatistics "${logSummary}"
    if [[ ${arrDetectorStatus[${detector}]} != "OK" ]]; then
      executePlanB
    fi
    echo "${detector} = ${arrDetectorStatus[${detector}]}" >> ${logFileShort}
  done

  #one email to the responsible - every time with a short summary
  if [[ -n ${MAILTO} && -n ${MAILshortSummary} ]]; then
    echo "mailing short summary to ${MAILTO}"
    mail -s "QA ready" ${MAILTO} < ${logFileShort}
  fi


  #remove lock
  rm -f ${lockFile}
  return 0
}

executePlanB()
{
  #in case of emergency
  #first check if we have the email of the detector expert defined,
  #if yes, append to the mailing list
  #NEEDS: $detector, $logSummary, $MAILTO, $MAILTO_DET

  local mailTo=${MAILTO}
  local detExpertEmailVar="MAILTO_${detector}"
  [[ -n "${!detExpertEmailVar}" ]] && mailTo+=" ${!detExpertEmailVar}"
  [[ -z ${mailTo} ]] && return 1

  echo
  echo "trouble detected, sending email to ${mailTo}"
  local mailoptions=""
  local file=""

  #attach the log summary
  file="${logSummary}"
  [[ ${MAILcompressLogs} == 1 && -f ${file} ]] && tar czf ${file}.tgz ${file} && file+=".tgz"
  [[ -f ${file} ]] && mailoptions+=" -a ${file}"

  #attach the crash plots
  file="${logDirectory}/stacktrace-log-${dateString}.png"
  [[ -f ${file} ]] && mailoptions+=" -a ${file}"
  file="${logDirectory}/stacktrace-core-${dateString}.png"
  [[ -f ${file} ]] && mailoptions+=" -a ${file}"

  #attach the full debug info
  if [[ -n "${MAILdebugInfo}" ]]; then
    file="${logDirectory}/stacktrace-log-${dateString}.tree"
    [[ ${MAILcompressLogs} == 1 && -f ${file} ]] && tar czf ${file}.tgz ${file} && file+=".tgz"
    [[ -f ${file} ]] && mailoptions+=" -a ${file}"
    file="${logDirectory}/stacktrace-core-${dateString}.tree"
    [[ ${MAILcompressLogs} == 1 && -f ${file} ]] && tar czf ${file}.tgz ${file} && file+=".tgz"
    [[ -f ${file} ]] && mailoptions+=" -a ${file}"
  fi

  #attach the full production log
  if [[ -n "${MAILfullProductionLog}" ]]; then
    file="${logFile}"
    [[ ${MAILcompressLogs} == 1 && -f ${file} ]] && tar czf ${file}.tgz ${file} && file+=".tgz"
    [[ -f ${file} ]] && mailoptions+=" -a ${file}"
  fi
  [[ -n ${MAILFROM} ]] && mailoptions+=" -r ${MAILFROM}"

  printLogStatistics ${logSummary} | mail -s "${detector} QA in need of assistance" ${mailoptions} ${mailTo}

  return 0
}

validate()
{
  summarizeLogs ${1}/* >> ${logSummary}
  logStatus=$?
  if [[ ${logStatus} -ne 0 ]]; then
    echo "WARNING not validated: ${1}"
    detectorStatus="planB"
    return 1
  fi
  return 0
}

substituteDetectorName()
{
  local det=$1
  local dir=$2
  [[ ${dir} =~ \%det ]] && det=${det,,} && echo ${dir/\%det/${det}}
  [[ ${dir} =~ \%DET ]] && det=${det} && echo ${dir/\%DET/${det}}
  return 0
}

main "$@"
