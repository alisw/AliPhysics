#!/bin/bash
# process QA output into plots and trending
# run without arguments for examples
# origin: Mikolaj Krzewicki, mkrzewic@cern.ch
#
if [ ${BASH_VERSINFO} -lt 4 ]; then
  echo "bash version >= 4 needed, you have ${BASH_VERSION}, exiting..."
  exit 1
fi

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

  #logging
  mkdir -p $logDirectory
  [[ ! -d $logDirectory ]] && echo "no log dir $logDirectory" && return 1
  logFile="$logDirectory/${0##*/}.${dateString}.log"
  touch ${logFile}
  [[ ! -f ${logFile} ]] && echo "cannot write logfile $logfile" && return 1
  echo "logFile = $logFile"

  #check lock
  lockFile=${workingDirectory}/runQA.lock
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
    #by default we expect the container in the QA root file to de named like
    #the detector
    detectorQAcontainerName=${detector}
    
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
    cd ${tmpDetectorRunDir}

    tmpPrefix=${tmpDetectorRunDir}/${outputDir}
    echo "running QA for ${detector}"
    echo "  outputDir=$outputDir"
    echo "  tmpPrefix=$tmpPrefix"
    
    #source the detector script
    #unset the detector functions from previous iterations (detectors)
    unset -f runLevelQA
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
      
      #first check if input file exists
      [[ ! -f ${inputFile%\#*} ]] && echo "file ${inputFile%\#*} not accessible" && continue

      if ! guessRunData ${inputFile}; then
        echo "could not guess run data from ${inputFile}"
        continue
      fi
      echo "anchorYear for ${originalPeriod} is: ${anchorYear}"

      tmpProductionDir=${tmpPrefix}/${dataType}/${year}/${period}/${pass}
      tmpRunDir=${tmpProductionDir}/000${runNumber}
      mkdir -p ${tmpRunDir}
      cd ${tmpRunDir}

      #check what kind of input file we have, default is a zip archive
      #set the inputs accordingly
      qaFile=""
      qaFileOuter=""
      highPtTree=""
      eventStatFile=""
      #it is possible we get the highPt trees from somewhere else
      #search the list of high pt trees for the proper run number
      if [[ -n ${inputListHighPtTrees} ]]; then
        highPtTree=$(egrep -m1 ${runNumber} ${inputListHighPtTrees})
        echo "loaded the highPtTree ${highPtTree} from external file ${inputListHighPtTrees}"
      fi
      #if we are explicit about the input file this takes precedence 
      #over earlier additions
      [[ "${inputFile}" =~ QAresults.root$ ]] && qaFile=${inputFile}
      [[ "${inputFile}" =~ QAresults_outer.root$ ]] && qaFileOuter=${inputFile}
      [[ "${inputFile}" =~ FilterEvents_Trees.root$ ]] && highPtTree=${inputFile}
      [[ "${inputFile}" =~ event_stat.root$ ]] && eventStatFile=${inputFile}
      if [[ "${inputFile}" =~ \.zip$ ]]; then
        [[ -z ${qaFile} ]] && qaFile=${inputFile}
        [[ -z ${qaFileOuter} ]] && qaFileOuter=${inputFile}
        [[ -z ${highPtTree} ]] && highPtTree=${inputFile}
        [[ -z ${eventStatFile} ]] && eventStatFile=${inputFile}
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
     
      echo qaFile=$qaFile
      echo qaFileOuter=$qaFileOuter
      echo highPtTree=$highPtTree
      echo eventStatFile=$eventStatFile
      echo ocdbStorage=${ocdbStorage}
      echo

      #standard QA based on QAresults.root file (and variants)
      if [[ -n ${qaFile} && $(type -t runLevelQA) =~ "function" ]]; then
        echo running ${detector} runLevelQA for run ${runNumber} from ${qaFile}
        ( runLevelQA "${qaFile}" ) &>> runLevelQA.log
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

      #perform some default actions:
      #if trending.root not created, create a default one
      if [[ ! -f trending.root ]]; then
        aliroot -b -q -l "$ALICE_ROOT/PWGPP/macros/simpleTrending.C(\"${qaFile}\",${runNumber},\"${detectorQAcontainerName}\",\"trending.root\",\"trending\",\"recreate\")" 2>&1 | tee -a runLevelQA.log
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
      if /bin/ls 000*/trending.root &>/dev/null; then
        hadd trending.root 000*/trending.root &> periodLevelQA.log
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
        echo "${HOSTNAME} ${dateString}" > ${periodLevelLock}
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
        planB=1
      fi

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
  return 0
}

executePlanB()
{
  #in case of emergency
  #first check if we have the email of the detector expert defined,
  #if yes, append to the mailing list
  local mailTo=${MAILTO}
  local detExpertEmailVar="MAILTO_${detector}"
  [[ -n "${!detExpertEmailVar}" ]] && mailTo+=" ${!detExpertEmailVar}"
  if [[ -n ${mailTo} ]]; then 
    echo
    echo "trouble detected, sending email to ${mailTo}"
    cat ${logSummary} | mail -s "${detector} QA in need of assistance" ${mailTo}
  fi
  return 0
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
    [[ ! -f ${log} ]] && continue
    errorSummary=$(validateLog ${log})
    validationStatus=$?
    [[ validationStatus -ne 0 ]] && logstatus=1
    if [[ ${validationStatus} -eq 0 ]]; then 
      #in pretend mode randomly report an error in rec.log some cases
      if [[ -n ${pretend} && "${log}" == "rec.log" ]]; then
        [[ $(( ${RANDOM}%2 )) -ge 1 ]] && echo "${log} BAD random error" || echo "${log} OK"
      else
        echo "${log} OK"
      fi
    elif [[ ${validationStatus} -eq 1 ]]; then
      echo "${log} BAD ${errorSummary}"
    elif [[ ${validationStatus} -eq 2 ]]; then
      echo "${log} OK MWAH ${errorSummary}"
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
  runMap="
  2010 108350 139517
  2011 140441 170593
  2012 171590 193766
  2013 194482 197692
  "

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
    echo "${var}=${value}"
    export ${var}="${value}"
  done
  return 0
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
  originalPass=""
  originalPeriod=""
  anchorYear=""

  shortRunNumber=""
  oldIFS=${IFS}
  local IFS="/"
  declare -a path=( $1 )
  IFS="${oldIFS}"
  local dirDepth=$(( ${#path[*]}-1 ))
  i=0
  for ((x=${dirDepth};x>=0;x--)); do

    [[ $((x-1)) -ge 0 ]] && local fieldPrev=${path[$((x-1))]}
    local field=${path[${x}]}
    local fieldNext=${path[$((x+1))]}

    [[ ${field} =~ ^[0-9]*$ && ${fieldNext} =~ (.*\.zip$|.*\.root$) ]] && legoTrainRunNumber=${field}
    [[ -n ${legoTrainRunNumber} && -z ${pass} ]] && pass=${fieldPrev}
    [[ ${field} =~ ^LHC[0-9][0-9][a-z].*$ ]] && period=${field%_*} && originalPeriod=${field}
    [[ ${field} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && runNumber=${field#000}
    [[ ${field} =~ ^[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && shortRunNumber=${field}
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field}
    [[ ${field} =~ ^(^sim$|^data$) ]] && dataType=${field}
    (( i++ ))
  done
  originalPass=${pass}
  [[ -n ${shortRunNumber} && "${legoTrainRunNumber}" =~ ${shortRunNumber} ]] && legoTrainRunNumber=""
  [[ -z ${legoTrainRunNumber} ]] && pass=${path[$((dirDepth-1))]}
  [[ "${dataType}" =~ ^sim$ ]] && pass="passMC" && runNumber=${shortRunNumber} && originalPass="" #for MC not from lego, the runnumber is identified as lego train number, thus needs to be nulled
  [[ -n ${legoTrainRunNumber} ]] && pass+="_lego${legoTrainRunNumber}"
  
  #modify the OCDB: set the year
  if [[ ${dataType} =~ sim ]]; then 
    anchorYear=$(run2year $runNumber)
    if [[ -z "${anchorYear}" ]]; then
      echo "WARNING: anchorYear not available for this production: ${originalPeriod}, runNumber: ${runNumber}. Cannot set the OCDB."
      return 1
    fi
    ocdbStorage=$(setYear ${anchorYear} ${ocdbStorage})
  else
    ocdbStorage=$(setYear ${year} ${ocdbStorage})
  fi

  #if [[ -z ${dataType} || -z ${year} || -z ${period} || -z ${runNumber}} || -z ${pass} ]];
  if [[ -z ${runNumber} ]]
  then
    #error condition
    return 1
  fi
  
  #ALL OK
  return 0
}

run2year()
{
  #for a given run print the year.
  #the run-year table is ${runMap} (a string)
  #defined in the config file
  #one line per year, format: year runMin runMax
  local run=$1
  [[ -z ${run} ]] && return 1
  local year=""
  local runMin=""
  local runMax=""
  while read year runMin runMax; do
    [[ -z ${year} || -z ${runMin} || -z ${runMax} ]] && continue
    [[ ${run} -ge ${runMin} && ${run} -le ${runMax} ]] && echo ${year} && break
  done < <(echo "${runMap}")
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

setYear()
{
  #set the year
  #  ${1} - year to be set
  #  ${2} - where to set the year
  local year1=$(guessYear ${1})
  local year2=$(guessYear ${2})
  local path=${2}
  [[ ${year1} -ne ${year2} && -n ${year2} && -n ${year1} ]] && path=${2/\/${year2}\//\/${year1}\/}
  echo ${path}
  return 0
}

guessYear()
{
  #guess the year from the path, pick the rightmost one
  local IFS="/"
  declare -a pathArray=( ${1} )
  local field
  local year
  for field in ${pathArray[@]}; do
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field}
  done
  echo ${year}
  return 0
}

hostInfo(){
#
# Hallo world -  Print AliRoot/Root/Alien system info
#

#
# HOST info
#
    echo --------------------------------------
        echo 
        echo HOSTINFO
        echo 
        echo HOSTINFO HOSTNAME"      "$HOSTNAME
        echo HOSTINFO DATE"          "`date`
        echo HOSTINFO gccpath"       "`which gcc` 
        echo HOSTINFO gcc version"   "`gcc --version | grep gcc`
        echo --------------------------------------    

#
# ROOT info
#
        echo --------------------------------------
        echo
        echo ROOTINFO
        echo 
        echo ROOTINFO ROOT"           "`which root`
        echo ROOTINFO VERSION"        "`root-config --version`
        echo 
        echo --------------------------------------


#
# ALIROOT info
#
        echo --------------------------------------
        echo
        echo ALIROOTINFO
        echo 
        echo ALIROOTINFO ALIROOT"        "`which aliroot`
        echo ALIROOTINFO VERSION"        "`echo $ALICE_LEVEL`
        echo ALIROOTINFO TARGET"         "`echo $ALICE_TARGET`
        echo 
        echo --------------------------------------

#
# Alien info
#
#echo --------------------------------------
#echo
#echo ALIENINFO
#for a in `alien --printenv`; do echo ALIENINFO $a; done 
#echo
#echo --------------------------------------

#
# Local Info
#
        echo PWD `pwd`
        echo Dir 
        ls -al
        echo
        echo
        echo
  
  return 0
}

main "$@"
