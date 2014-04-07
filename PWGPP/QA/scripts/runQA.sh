#!/bin/bash
main()
{
  if [[ -z $1 ]]; then
    echo "Usage: "
    echo "  ${0##*/} option=value [option=value]"
    echo "  at least inputList should be specified, or configFile containing it:"
    echo "  ${0##*/} inputList=file.list"
    echo "  options override config file (if any), e.g.:"
    echo "  ${0##*/} configFile=runQA.config inputList=file.list outputDirectory=%det"
    return 1
  fi
 
  if ! parseConfig $@; then
    ${0}
    return 1
  fi

  [[ -z $ALICE_ROOT ]] && echo "ALICE_ROOT not defined" && return 1

  ocdbregex='raw://'
  if [[ ${ocdbStorage} =~ ${ocdbregex} ]]; then
    alien-token-init ${alienUserName}
    #this is a hack! alien-token init seems not enough
    #but the gclient_env script messes up the LD_LIBRARY_PATH
    while read x; do
      eval ${x};
    done < <(grep -v "LD_LIBRARY_PATH" /tmp/gclient_env_${UID})
  fi

  updateQA $@
}

updateQA()
{
  umask 0002
  parseConfig $@

  #be paranoid and make some full paths
  [[ ! -f ${inputList} ]] && echo "no input list: ${inputList}" && return 1
  inputList=$(get_realpath ${inputList})
  mkdir -p ${workingDirectory}
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

  dateString=$(date +%Y-%m-%d-%H-%M)
  echo "Start time QA process: $dateString"

  #logging
  mkdir -p $logDirectory
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && return 1
  logFile="$logDirectory/${0##*/}.${dateString}.log"
  touch ${logFile}
  [[ ! -f ${logFile} ]] && echo "cannot write logfile $logfile" && return 1
  echo "logFile = $logFile"

  #check lock
  lockFile=${logDirectory}/runQA.lock
  [[ -f ${lockFile} ]] && echo "lock ${lockFile} exists!" | tee ${logFile} && return 1
  touch ${lockFile}
  [[ ! -f ${lockFile} ]] && echo "cannot lock $lockFile" | tee ${logFile} && return 1
  
  exec &>${logFile}

  ################################################################
  #ze detector loop
  for detectorScript in $ALICE_ROOT/PWGPP/QA/detectorQAscripts/*; do
    echo
    echo "##############################################"
    echo $(date)
    unset planB
    [[ ! ${detectorScript} =~ .*\.sh$ ]] && continue
    detector=${detectorScript%.sh}
    detector=${detector##*/}
    
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
    outputDir=$(substituteDetectorName ${detector} ${outputDirectory})
    tmpDetectorRunDir=${workingDirectory}/tmpQAtmpRunDir${detector}-${dateString}
    if ! mkdir -p ${tmpDetectorRunDir}; then
      echo "cannot create the temp dir $tmpDetectorRunDir"
      continue
    fi
    cd ${tmpDetectorRunDir}

    tmpPrefix=${tmpDetectorRunDir}/${outputDir}
    echo "running QA for ${detector}"
    echo "  outputDir=$outputDir"
    echo "  tmpPrefix=$tmpPrefix"
    
    #unset the detector functions from previous iterations (detectors)
    unset -f runLevelQA
    unset -f periodLevelQA
    unset -f runLevelHighPtTreeQA
    unset -f periodLevelHighPtTreeQA
    source ${detectorScript}

    #################################################################
    #produce the QA and trending tree for each file (run)
    unset arrOfTouchedProductions
    declare -A arrOfTouchedProductions
    while read qaFile; do
      echo
      echo $(date)
      
      #first check if input file exists
      [[ ! -f ${qaFile%\#*} ]] && echo "file ${qaFile%\#*} not accessible" && continue

      if ! guessRunData ${qaFile}; then
        echo "could not guess run data from ${qaFile}"
        continue
      fi

      tmpProductionDir=${tmpPrefix}/${dataType}/${year}/${period}/${pass}
      tmpRunDir=${tmpProductionDir}/000${runNumber}
      mkdir -p ${tmpRunDir}
      cd ${tmpRunDir}

      #by default we expect to have everything in the same archive
      highPtTree=${qaFile}

      #maybe the input is not an archive, but a file
      [[ "${qaFile}" =~ QAresults.root$ ]] && highPtTree=""
      [[ "${qaFile}" =~ FilterEvents_Trees.root$ ]] && qaFile=""

      #it is possible we get the highPt trees from somewhere else
      #search the list of high pt trees for the proper run number
      if [[ -n ${inputListHighPtTrees} ]]; then
        highPtTree=$(egrep -m1 ${runNumber} ${inputListHighPtTrees})
        echo "loaded the highPtTree ${highPtTree} from external file ${inputListHighPtTrees}"
      fi
      
      echo qaFile=$qaFile
      echo highPtTree=$highPtTree

      #what if we have a zip archive?
      if [[ "$qaFile" =~ .*.zip$ ]]; then
        if unzip -l ${qaFile} | egrep "QAresults.root" &>/dev/null; then
          qaFile="${qaFile}#QAresults.root"
        else
          qaFile=""
        fi
      fi
      if [[ "$highPtTree" =~ .*.zip$ ]]; then
        if unzip -l ${highPtTree} | egrep "FilterEvents_Trees.root" &>/dev/null; then
          highPtTree="${highPtTree}#FilterEvents_Trees.root"
        else
          highPtTree=""
        fi
      fi
     
      if [[ -n ${qaFile} && $(type -t runLevelQA) =~ "function" ]]; then
        echo running ${detector} runLevelQA for run ${runNumber} from ${qaFile}
        runLevelQA "${qaFile}" &> runLevelQA.log
        #perform some default actions:
        #if trending.root not created, create a default one
        if [[ ! -f trending.root ]]; then
          aliroot -b -q -l "$ALICE_ROOT/PWGPP/macros/simpleTrending.C(\"${qaFile}\",${runNumber},\"${detector}\",\"trending.root\",\"trending\",\"recreate\")" 2>&1 | tee -a runLevelQA.log
        fi
        if [[ -f trending.root ]]; then
          arrOfTouchedProductions[${tmpProductionDir}]=1
        else
          echo "trending.root not created"
        fi
      fi
      #expert QA based on high pt trees
      if [[ -n ${highPtTree} && $(type -t runLevelHighPtTreeQA) =~ "function" ]]; then
        echo running ${detector} runLevelHighPtTreeQA for run ${runNumber} from ${highPtTree}
        runLevelHighPtTreeQA "${highPtTree}" &> runLevelHighPtTreeQA.log
        arrOfTouchedProductions[${tmpProductionDir}]=1
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
        if ! guessRunData "${dir}/dummyName"; then
          echo "could not guess run data from ${dir}"
          continue
        fi

        #before moving - VALIDATE!!!
        if ! validate ${dir}; then 
          continue
        fi

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
      echo tmpPeriodLevelQAdir="${tmpProductionDir}/periodLevelQA"
      if ! mkdir -p ${tmpPeriodLevelQAdir}; then continue; fi
      cd ${tmpPeriodLevelQAdir}

      #link the final list of per-run dirs here, just the dirs
      #to have a clean working directory
      unset linkedStuff
      declare -a linkedStuff
      for x in ${productionDir}/000*; do [[ -d $x ]] && ln -s $x && linkedStuff+=(${x##*/}); done

      #merge trending files if any
      if /bin/ls 000*/trending.root &>/dev/null; then
        hadd trending.root 000*/trending.root &> periodLevelQA.log
      fi
      
      #run the period level trending/QA
      if [[ -f "trending.root" && $(type -t periodLevelQA) =~ "function" ]]; then
        echo running ${detector} periodLevelQA for production ${period}/${pass}
        periodLevelQA trending.root &>> periodLevelQA.log
      else 
        echo "WARNING: not running ${detector} periodLevelQA for production ${period}/${pass}, no trending.root"
      fi

      if ! validate ${PWD}; then continue; fi

      #here we are validated so move the produced QA to the final place
      #clean up linked stuff first
      [[ -n ${linkedStuff[@]} ]] && rm ${linkedStuff[@]}
      #some of the output could be a directory, so handle that
      #TODO: maybe use rsync?
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

      #remove the temp dir
      rm -rf ${tmpPeriodLevelQAdir}
    
    done

    cd ${workingDirectory}

    if [[ -z ${planB} ]]; then
      echo
      echo removing ${tmpDetectorRunDir}
      rm -rf ${tmpDetectorRunDir}
    else
      executePlanB
    fi
  done #end of detector loop

  #remove lock
  rm -f ${lockFile}
}

executePlanB()
{
  #in case of emergency
  if [[ -n ${MAILTO} ]]; then 
    echo
    echo "trouble detected, sending email to ${MAILTO}"

    grep BAD ${logSummary} | mail -s "qa in need of assistance" ${MAILTO}
  fi
}

validate()
{
  summarizeLogs ${1} >> ${logSummary}
  logStatus=$?
  if [[ ${logStatus} -ne 0 ]]; then 
    echo "WARNING not validated: ${1}"
    planB=1
    return 1
  fi
  return 0
}

summarizeLogs()
{
  local dir=$1
  [[ ! -d ${dir} ]] && dir=${PWD}

  #print a summary of logs
  logFiles=(
      "*.log"
      "stdout"
      "stderr"
  )

  #check logs
  local logstatus=0
  for log in ${dir}/${logFiles[*]}; do
    finallog=${PWD%/}/${log}
    [[ ! -f ${log} ]] && continue
    errorSummary=$(validateLog ${log})
    validationStatus=$?
    [[ validationStatus -ne 0 ]] && logstatus=1
    if [[ ${validationStatus} -eq 0 ]]; then 
      #in pretend mode randomly report an error in rec.log some cases
      if [[ -n ${pretend} && "${log}" == "rec.log" ]]; then
        [[ $(( ${RANDOM}%2 )) -ge 1 ]] && echo "${finallog} BAD random error" || echo "${finallog} OK"
      else
        echo "${finallog} OK"
      fi
    elif [[ ${validationStatus} -eq 1 ]]; then
      echo "${finallog} BAD ${errorSummary}"
    elif [[ ${validationStatus} -eq 2 ]]; then
      echo "${finallog} OK MWAH ${errorSummary}"
    fi
  done

  #report core files
  while read x; do
    echo ${x}
    chmod 644 ${x}
    gdb --batch --quiet -ex "bt" -ex "quit" aliroot ${x} > stacktrace_${x//\//_}.log
  done < <(/bin/ls ${PWD}/*/core 2>/dev/null; /bin/ls ${PWD}/core 2>/dev/null)

  return ${logstatus}
}

validateLog()
{
  log=${1}
  errorConditions=(
            'There was a crash'
            'floating'
            'error while loading shared libraries'
            'std::bad_alloc'
            's_err_syswatch_'
            'Thread [0-9]* (Thread'
            'AliFatal'
            'core dumped'
            '\.C.*error:.*\.h: No such file'
            'segmentation'
            'Interpreter error recovered'
  )

  warningConditions=(
            'This is serious'
  )

  local logstatus=0
  local errorSummary=""
  local warningSummary=""

  for ((i=0; i<${#errorConditions[@]};i++)); do
    local tmp=$(grep -m1 -e "${errorConditions[${i}]}" ${log})
    [[ -n ${tmp} ]] && tmp+=" : "
    errorSummary+=${tmp}
  done

  for ((i=0; i<${#warningConditions[@]};i++)); do
    local tmp=$(grep -m1 -e "${warningConditions[${i}]}" ${log})
    [[ -n ${tmp} ]] && tmp+=" : "
    warningSummary+=${tmp}
  done

  if [[ -n ${errorSummary} ]]; then 
    echo "${errorSummary}"
    return 1
  fi

  if [[ -n ${warningSummary} ]]; then
    echo "${warningSummary}"
    return 2
  fi

  return 0
}

parseConfig()
{
  args=("$@")

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
  #email to
  #MAILTO="fbellini@cern.ch"

  #first, check if the config file is configured
  #is yes - source it so that other options can override it
  #if any
  for opt in "${args[@]}"; do
    if [[ ${opt} =~ configFile=.* ]]; then
      eval "${opt}"
      [[ ! -f ${configFile} ]] && echo "configFile ${configFile} not found, exiting..." && return 1
      echo "using config file: ${configFile}"
      source "${configFile}"
      break
    fi
  done

  #then, parse the options as they override the options from file
  for opt in "${args[@]}"; do
    if [[ ! "${opt}" =~ .*=.* ]]; then
      echo "badly formatted option ${var}, should be: option=value, stopping..."
      return 1
    fi
    local var="${opt%%=*}"
    local value="${opt#*=}"
    echo "${var} = ${value}"
    export ${var}="${value}"
  done
}

guessRunData()
{
  #guess the period from the path, pick the rightmost one
  period=""
  runNumber=""
  year=""
  pass=""
  legoTrainRunNumber=""
  dataType=""

  local shortRunNumber=""
  local IFS="/"
  declare -a path=( $1 )
  local dirDepth=$(( ${#path[*]}-1 ))
  i=0
  for ((x=${dirDepth};x>=0;x--)); do

    [[ $((x-1)) -ge 0 ]] && local fieldPrev=${path[$((x-1))]}
    local field=${path[${x}]}
    local fieldNext=${path[$((x+1))]}

    [[ ${field} =~ ^[0-9]*$ && ${fieldNext} =~ (.*\.zip$|.*\.root$) ]] && legoTrainRunNumber=${field}
    [[ -n ${legoTrainRunNumber} && -z ${pass} ]] && pass=${fieldPrev}
    [[ ${field} =~ ^LHC[0-9][0-9][a-z].*$ ]] && period=${field%_*}
    [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && runNumber=${field#000}
    [[ ${field} =~ ^[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && shortRunNumber=${field}
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field}
    [[ ${field} =~ ^(^sim$|^data$) ]] && dataType=${field}
    (( i++ ))
  done
  [[ -z ${legoTrainRunNumber} ]] && pass=${path[$((dirDepth-1))]}
  [[ "${dataType}" =~ ^sim$ ]] && pass="passMC" && runNumber=${shortRunNumber}
  
  #if [[ -z ${dataType} || -z ${year} || -z ${period} || -z ${runNumber}} || -z ${pass} ]];
  if [[ -z ${runNumber}} ]];
  then
    #error condition
    return 1
  else
    #ALL OK
    return 0
  fi
}

substituteDetectorName()
{
  local det=$1
  local dir=$2
  [[ ${dir} =~ \%det ]] && det=${det,,} && echo ${dir/\%det/${det}}
  [[ ${dir} =~ \%DET ]] && det=${det} && echo ${dir/\%DET/${det}}
}

get_realpath() 
{
  if [[ -f "$1" ]]
  then
    # file *must* exist
    if cd "$(echo "${1%/*}")" &>/dev/null
    then
      # file *may* not be local
      # exception is ./file.ext
      # try 'cd .; cd -;' *works!*
      local tmppwd="$PWD"
      cd - &>/dev/null
    else
      # file *must* be local
      local tmppwd="$PWD"
    fi
  else
    # file *cannot* exist
    return 1 # failure
  fi
  # reassemble realpath
  echo "$tmppwd"/"${1##*/}"
  return 0 # success
}

main $@
