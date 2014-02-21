#!/bin/bash
main()
{
  if [[ -z $1 ]]; then
    echo "Usage: $0 configFile [optionalStuff]"
    echo "  - optionalStuff overrides config file, e.g.:"
    echo "       $0 configFile inputList=somefile.list outputDirectory='${PWD}'/output"
    exit 1
  fi
  [[ ! -f $1 ]] && echo "argument not a file" && exit 1  

 
  configFile=${1}
  shift 1

  [[ -z $ALICE_ROOT ]] && source $alirootEnv
  [[ -z $ALICE_ROOT ]] && echo "ALICE_ROOT not defined" && exit 1

  ocdbregex='raw://'
  [[ ${ocdbStorage} =~ ${ocdbregex} ]] && alien-token-init

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
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && exit 1
  logFile="$logDirectory/${0##*/}.${dateString}.log"
  touch ${logFile}
  [[ ! -f ${logFile} ]] && echo "cannot write logfile $logfile" && exit 1
  echo "logFile = $logFile"
  exec &>${logFile}

  #check lock
  lockFile=${logDirectory}/runQA.lock
  [[ -f ${lockFile} ]] && echo "lock ${lockFile} exists!" && exit 1
  touch ${lockFile}
  [[ ! -f ${lockFile} ]] && echo "cannot lock $lockFile" && exit 1
  
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
    outputDir=$(substituteDetectorName ${detector} ${outputDirectory})
    tmpRunDir=${workingDirectory}/tmpQArunDir${detector}
    if ! mkdir -p ${tmpRunDir}; then
      echo "cannot create the temp dir $tmpRunDir"
      continue
    fi
    cd ${tmpRunDir}

    tmpPrefix=${tmpRunDir}/${outputDir}
    echo "running QA for ${detector}"
    echo "  outputDir=$outputDir"
    echo "  tmpPrefix=$tmpPrefix"
    
    unset -f runLevelQA
    unset -f periodLevelQA
    source ${detectorScript}

    #################################################################
    #produce the QA and trending tree for each file (run)
    while read qaFile; do
    
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
    
      productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}
      
      echo mkdir -p ${productionDir}
      mkdir -p ${productionDir}
      if [[ ! -d ${productionDir} ]]; then 
        echo "cannot make productionDir $productionDir" && continue
      fi
      
      #move to final destination
      for dir in ${tmpProductionDir}/*; do
        oldRunDir=${outputDir}/${dir#${tmpPrefix}}
        if [[ -d ${oldRunDir} ]]; then
          echo "removing old run: rm -rf ${oldRunDir}"
          rm -rf ${oldRunDir}
        fi
        echo "moving output to final destination"
        echo mv -f ${dir} ${productionDir}
        mv -f ${dir} ${productionDir}
      done
    
      guessRunData "${productionDir}/dummyName"

      cd ${productionDir}

      rm -f trending.root
      hadd trending.root 000*/trending.root

      echo running ${detector} periodLevelQA for production ${period}/${pass}
      periodLevelQA trending.root &> periodLevelQA.log
      
      cd ${tmpRunDir}
    
    done

    cd ${workingDirectory}
    echo cleaning up: rm -rf ${tmpRunDir}
    rm -rf ${tmpRunDir}
  done

  #remove lock
  rm -f ${lockFile}
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
