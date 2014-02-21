#!/bin/bash
main()
{
  if [[ -z $1 ]]; then
    echo "Usage: $0 configFile [optionalStuff]"
    echo "  - optionalStuff overrides config file, e.g.:"
    echo "       $0 configFile inputList=somefile.list outputDirectory='${PWD}'/output"
    return 1
  fi
  [[ ! -f $1 ]] && echo "argument not a file" && return 1  
 
  configFile=${1}
  shift 1
  parseConfig ${configFile} $@

  [[ -z $ALICE_ROOT ]] && source ${alirootEnv}
  [[ -z $ALICE_ROOT ]] && echo "ALICE_ROOT not defined" && return 1

  ocdbregex='raw://'
  if [[ ${ocdbStorage} =~ ${ocdbregex} ]]; then
    
    alien-token-init
  fi

  updateQA ${configFile} $@
}

updateQA()
{
  #this guy takes config file as only positional argument
  #optional stugg allowerd to override config
  umask 0002
  configFile=${1}
  shift 1
  parseConfig $configFile $@

  dateString=$(date +%Y-%m-%d-%H-%M)
  echo "Start time QA process: $dateString"

  #logging
  mkdir -p $logDirectory
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && return 1
  logFile="$logDirectory/${0##*/}.${dateString}.log"
  touch ${logFile}
  [[ ! -f ${logFile} ]] && echo "cannot write logfile $logfile" && return 1
  echo "logFile = $logFile"
  exec &>${logFile}

  #check lock
  lockFile=${logDirectory}/runQA.lock
  [[ -f ${lockFile} ]] && echo "lock ${lockFile} exists!" && return 1
  touch ${lockFile}
  [[ ! -f ${lockFile} ]] && echo "cannot lock $lockFile" && return 1
  
  #be paranoid and make some full paths
  inputList=$(readlink -f ${inputList})
  workingDirectory=$(readlink -f ${workingDirectory})
  mkdir -p ${workingDirectory}
  if [[ ! -d ${workingDirectory} ]]; then
    echo "working dir $workingDirectory does not exist and cannot be created"
    return 1
  fi
  cd ${workingDirectory}

  ################################################################
  #ze detector loop
  for detectorScript in $ALICE_ROOT/PWGPP/QA/detectorQAscripts/*; do

    [[ ! ${detectorScript} =~ .*\.sh ]] && continue
    detector=${detectorScript%.sh}
    detector=${detector##*/}
    
    #skip if excluded
    skipDetector=0
    for excluded in ${excludeDetectors}; do
      if [[ ${detector} =~ ${excluded} ]]; then
        echo "${detector} is excluded in config, skipping..."
        skipDetector=1
        break
      fi
    done
    [[ ${skipDetector} -eq 1 ]] && continue

    logSummary=${logDirectory}/summary-${detector}-${dateString}.log
    outputDir=$(substituteDetectorName ${detector} ${outputDirectory})
    tmpRunDir=${workingDirectory}/tmpQArunDir${detector}
    if ! mkdir -p ${tmpRunDir}; then
      echo "cannot create the temp dir $tmpRunDir"
      continue
    fi
    cd ${tmpRunDir}

    tmpPrefix=${tmpRunDir}/${outputDir}
    echo
    echo "##############################################"
    echo "running QA for ${detector}"
    echo "  outputDir=$outputDir"
    echo "  tmpPrefix=$tmpPrefix"
    
    unset -f runLevelQA
    unset -f periodLevelQA
    source ${detectorScript}

    #################################################################
    #produce the QA and trending tree for each file (run)
    while read qaFile; do
      echo

      guessRunData ${qaFile}

      runDir=${tmpPrefix}/${dataType}/${year}/${period}/${pass}/000${runNumber}
      mkdir -p ${runDir}
      cd ${runDir}

      #handle the case of a zip archive
      [[ "$qaFile" =~ .*.zip$ ]] && qaFile="${qaFile}#QAresults.root"
      
      echo running ${detector} runLevelQA for run ${runNumber} from ${qaFile}
      runLevelQA ${qaFile} &> runLevelQA.log

      cd ${tmpRunDir}
    
    done < ${inputList}

    #################################################################
    #cache which productions were (re)done
    arrOfTouchedProductions=( $(/bin/ls -d ${tmpPrefix}/*/*/*/*) )
    echo "list of processed productions:"
    echo "    ${arrOfTouchedProductions[@]}"
    echo
    #################################################################
    #(re)do the merging/trending in the final destination
    for tmpProductionDir in ${arrOfTouchedProductions[@]}; do
      echo
      echo "running period level stuff in ${tmpProductionDir}"
    
      productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}

      mkdir -p ${productionDir}
      if [[ ! -d ${productionDir} ]]; then 
        echo "cannot make productionDir $productionDir" && continue
      fi
      cd ${productionDir}
      
      #move to final destination
      for dir in ${tmpProductionDir}/*; do
        oldRunDir=${outputDir}/${dir#${tmpPrefix}}
        guessRunData "${dir}/dummyName"

        #before moving - VALIDATE!!!
        if ! validate ${dir}; then continue; fi

        #summarizeLogs ${dir} >> ${logSummary}
        #logStatus=$?
        #if [[ ${logStatus} -ne 0 ]]; then 
        #  echo "WARNING: run not validated: ${dir}"
        #  planB=1
        #  continue
        #fi

        if [[ -d ${oldRunDir} ]]; then
          echo "removing old ${period}/${pass}/${runNumber}"
          rm -rf ${oldRunDir}
        fi
        echo "moving new ${runNumber} to ${productionDir}"
        mv -f ${dir} ${productionDir}
      done
    
      echo running ${detector} periodLevelQA for production ${period}/${pass}
      rm -f trending.root
      hadd trending.root 000*/trending.root &> periodLevelQA.log
      periodLevelQA trending.root &>> periodLevelQA.log
      
      if ! validate ${PWD}; then continue; fi
      #summarizeLogs ${PWD} >> ${logSummary}
      #logStatus=$?
      #if [[ ${logStatus} -ne 0 ]]; then 
      #  echo "WARNING: period ${period} not validated: ${dir}"
      #  planB=1
      #  continue
      #fi

      cd ${tmpRunDir}
    
    done

    cd ${workingDirectory}

    if [[ -z ${planB} ]]; then
      echo
      echo removing ${tmpRunDir}
      rm -rf ${tmpRunDir}
    else
      executePlanB
    fi
  done

  #remove lock
  rm -f ${lockFile}
}

executePlanB()
{
  #in case of emergency
  if [[ -n ${MAILTO} ]]; then 
    echo
    echo "trouble detected, sending email to ${MAILTO}"

    cat ${logSummary} | mail -s "qa in need of assistance" ${MAILTO}
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
  [[ -z ${dir} ]] && dir=${PWD}

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
  configFile=${1}
  shift

  #first, source the config file
  if [ -f ${configFile} ]; then
    source ${configFile}
  else
    echo "config file ${configFile} not found!, skipping..."
  fi

  #then, parse the options as they override the options from file
  while [[ -n ${1} ]]; do
    local var=${1#--}
    eval "${var}"
    shift
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
  
  if [[ -z ${dataType} || -z ${year} || -z ${period} || -z ${runNumber}} || -z ${pass} ]];
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

main $@
