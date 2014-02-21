#!/bin/bash
main()
{
  runQA $@
}

runQA()
{
  umask 0002
  dateString=$(date +%Y-%m-%d-%H-%M)
  [[ -z $1 ]] && echo "Usage: $0 configFile" && exit 1
  [[ ! -f $1 ]] && echo "argument not a file" && exit 1
  
  configFile=${1}
  shift 1
  parseConfig $configFile $@

  #be paranoid and make workingDirectory a full path
  workingDirectory=$(readlink -f ${workingDirectory})
  cd ${workingDirectory}

  [[ -z $ALICE_ROOT ]] && source $alirootEnv
  [[ -z $ALICE_ROOT ]] && echo "ALICE_ROOT not defined" && exit 1

  mkdir -p $logDirectory
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && exit 1

  logFile="$logDirectory/${0##*/}.$dateString.log"
  touch $logFile
  [[ ! -f $logFile ]] && echo "cannot write logfile $logfile" && exit 1
  exec 1>$logFile 2>&1

  dateString=$(date +%Y-%m-%d-%H-%M)
  echo "Start time QA process: $dateString"

  for detectorScript in $ALICE_ROOT/PWGPP/QA/detectorQAscripts/*; do

    [[ ! ${detectorScript} =~ .*\.sh ]] && continue
    detector=${detectorScript%.sh}
    detector=${detector##*/}
    outputDir=$(substituteDetectorName ${detector} ${outputDirectory})
    tmpPrefix=${workingDirectory}/runDir/${outputDir}
    echo outputDir=$outputDir
    echo tmpPrefix=$tmpPrefix
    echo detector=$detector
    
    unset -f runLevelQA
    unset -f periodLevelQA
    source ${detectorScript}

    #################################################################
    #produce the QA and trending tree for each file (run)
    while read qaFile; do
    
      guessRunData ${qaFile}

      runDir=${tmpPrefix}/${dataType}/${year}/${period}/${pass}/000${runNumber}
      echo runDir=$runDir
      mkdir -p ${runDir}
      cd ${runDir}

      [[ "$qaFile" =~ .*.zip$ ]] && qaFile="${qaFile}#QAresults.root"
      
      runLevelQA ${qaFile}

      cd ${workingDirectory}
    
    done < ${inputList}

    #################################################################
    #cache which productions were (re)done
    arrOfTouchedProductions=( $(/bin/ls -d ${tmpPrefix}/*/*/*/*) )
    echo arrOfTouchedProductions=${arrOfTouchedProductions[@]}
    
    #################################################################
    #(re)do the merging/trending in the final destination
    for tmpProductionDir in ${arrOfTouchedProductions[@]}; do
    
      echo productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}
      productionDir=${outputDir}/${tmpProductionDir#${tmpPrefix}}
      
      echo mkdir -p ${productionDir}
      mkdir -p ${productionDir}
      
      echo ln -s ${tmpProductionDir}/* ${productionDir}
      ln -s ${tmpProductionDir}/* ${productionDir}
    
      guessRunData "${tmpProductionDir}/dummyName"

      cd ${productionDir}

      hadd trending.root 000*/trending.root

      periodLevelQA trending.root
      
      cd ${workingDirectory}
    
    done

  done

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
}

substituteDetectorName()
{
  local det=$1
  local dir=$2
  [[ ${dir} =~ \%det ]] && det=${det,,} && echo ${dir/\%det/${det}}
  [[ ${dir} =~ \%DET ]] && det=${det} && echo ${dir/\%DET/${det}}
}

main $@
