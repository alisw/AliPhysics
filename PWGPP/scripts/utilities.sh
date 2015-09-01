#!/usr/bin/env bash
#library of useful PWGPP related bash functions
#it REQUIRES BASH 4 !!!!
#blame: Mikolaj Krzewicki, mkrzewic@cern.ch

if [ ${BASH_VERSINFO} -lt 4 ]; then
  echo "bash version >= 4 needed, you have ${BASH_VERSION}, exiting..."
  exit 1
fi

PWGPP_runMap="
2010 108350 139517
2011 140441 170593
2012 171590 193766
2013 194308 199146
2014 202369 206695
2015 208505 999999
2016 999999 999999
"

parseConfig()
{
  # parse command line arguments, they have to be in the form
  #  option=value
  # they are then set in the environment
  # additionally another variable named: parseConfig__ORIGINAL__${option}
  # is set to have a fallback.
  # The prefix can be changed by specifying the
  # option originalOptionPrefix="some_prefix". Set to "" to switch off.
  #
  # optionally a config file can be specified in the arguments:
  #  configFile=<someFile>
  # config file sets variables: option=value
  # command line options override config file
  #
  # recommended way of using (at the beginning of your funcion/script):
  #  if ! parseConfig "${@}"; then return; fi
  
  local args=("$@")
  local opt=""
  local originalOptionPrefix="parseConfig__ORIGINAL__"
  
  #first check if we will need to decode spaces
  local encodedSpaces=""
  for opt in "${args[@]}"; do
    [[ "${opt}" =~ encodedSpaces=.* ]] \
      && encodedSpaces=1
    [[ "${opt}" =~ originalOptionPrefix=.* ]] \
      && originalOptionPrefix="${opt#*=}"
  done

  #then look for a configFile (if any)
  for opt in "${args[@]}"; do
    if [[ ${opt} =~ configFile=.* ]]; then
      eval "${opt}"
      [[ ! -f ${configFile} ]] \
        && echo "configFile ${configFile} not found, exiting..." \
        && return 1
      echo "using config file: ${configFile}"
      source "${configFile}"
      break
    fi
  done

  #then, parse the options as they override the options from configFile
  for opt in "${args[@]}"; do
    [[ -n ${encodedSpaces} ]] && opt="$(decSpaces ${opt})"
    if [[ ! "${opt}" =~ .*=.* ]]; then
      echo "badly formatted option ${var}, should be: option=value, stopping..."
      return 1
    fi
    local var="${opt%%=*}"
    local value="${opt#*=}"
    #echo "${var}=${value}"
    export ${var}="${value}"
    [[ -n ${originalOptionPrefix} ]] && export ${originalOptionPrefix}${var}="${value}"
  done
  return 0
}

guessRunData()
{
  #guess the period from the path, pick the rightmost one
  #if $ocdbStorage is set it will be reset to the anchorYear (for MC)
  period=""
  runNumber=""
  year=""
  pass=""
  legoTrainRunNumber=""
  dataType=""
  dataAOD=""
  originalPass=""
  originalPeriod=""
  anchorYear=""
  shortRunNumber=""

  local oldIFS=${IFS}
  local IFS="/"
  [[ -z ${1} ]] && return 1
  declare -a path=( $1 )
  IFS="${oldIFS}"
  local dirDepth=$(( ${#path[*]}-1 ))
  for ((x=${dirDepth};x>=0;x--)); do

    [[ $((x-2)) -ge 0 ]] && local fieldPrevPrev=${path[$((x-2))]}
    [[ $((x-1)) -ge 0 ]] && local fieldPrev=${path[$((x-1))]}
    local field=${path[${x}]}
    local fieldNext=${path[$((x+1))]}
    local fieldNextNext=${path[$((x+2))]}

    [[ -z ${legoTrainRunNumber} && ${field} =~ ^[0-9]*$ && ${fieldNext} =~ (.*\.zip$|.*\.root$) ]] && legoTrainRunNumber=${field}
    [[ -z ${runNumber} && ${fieldPrev} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && runNumber=${fieldPrev#000}
    [[ -z ${pass} && -n ${runNumber} ]] && pass=${field}
    [[ ${field} =~ ^LHC[0-9][0-9][a-z].*$ ]] && period=${field%_*} && originalPeriod=${field}
    [[ ${field} =~ ^[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && shortRunNumber=${field}
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year=${field}
    [[ ${field} =~ ^(^sim$|^data$) ]] && dataType=${field}
    [[ "${field}" =~ "${pass}" && ${fieldNext} =~ AOD ]] && legoTrainRunNumber="" && dataAOD="AOD"

    [[ ${field} =~ ^LHC[0-9][0-9][a-z].*$ && ${fieldPrev} =~ ^20[0-9][0-9]$ && ! ${fieldNext} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ ]] && pass=${fieldNext}

    [[ ${field} =~ "ESDs" && ${fieldPrev} =~ ^000[0-9][0-9][0-9][0-9][0-9][0-9]$ && ${fieldPrevPrev} =~ ^LHC[0-9][0-9][a-z].*$ ]] && pass=${fieldNext}

  done
  originalPass=${pass}

  if [[ ${dataType} =~ sim ]]; then
    [[ -n ${shortRunNumber} && -z ${runNumber} ]] && runNumber=${shortRunNumber}
    pass="passMC"
    originalPass="" #for MC not from lego, the runNumber is identified as lego train number, thus needs to be nulled
    anchorYear=$(run2year $runNumber)
    if [[ -z "${anchorYear}" ]]; then
      echo "WARNING: anchorYear not available for this production: ${originalPeriod}, runNumber: ${runNumber}. Cannot set the OCDB."
      return 1
    fi
    #modify the OCDB: set the year
    ocdbStorage=$(setYear ${anchorYear} ${ocdbStorage})
  else
    ocdbStorage=$(setYear ${year} ${ocdbStorage})
  fi

  [[ -n ${shortRunNumber} && -z ${runNumber} && -z {dataType} ]] && runNumber=${shortRunNumber}
  [[ -n ${shortRunNumber} && "${legoTrainRunNumber}" =~ ${shortRunNumber} ]] && legoTrainRunNumber=""
  [[ -n ${legoTrainRunNumber} ]] && pass+="_lego${legoTrainRunNumber}"
  [[ -n ${dataAOD} ]] && pass="${pass}_AOD"

  #if [[ -z ${dataType} || -z ${year} || -z ${period} || -z ${runNumber}} || -z ${pass} ]];
  #if [[ -z ${runNumber} ]]
  #then
  #  #error condition
  #  return 1
  #fi
  
  #ALL OK
  return 0
}

guessRunNumber()
(
  #guess the run number from the path, pick the rightmost one
  if guessRunData "${1}"; then
    echo ${runNumber}
    return 0
  fi
  return 1
)

guessYear()
(
  #guess the year from the path, pick the rightmost one
  if guessRunData "${1}"; then
    echo ${year}
    return 0
  fi
  return 1
)

guessPeriod()
(
  #guess the period from the path, pick the rightmost one
  if guessRunData "${1}"; then
    echo ${period}
    return 0
  fi
  return 1
)

setYear()
{
  #set the year in the string
  #usualy used to modify the year in $ocdbStorage
  #  ${1} - year to be set
  #  ${2} - where to set the year
  #if AUTOYEAR is present in target - it will be replaced by the year
  local yearSource=$(guessYearFast ${1})
  local yearTarget=$(guessYearFast ${2})
  local path=${2}
  [[ ${yearSource} -ne ${yearTarget} && -n ${yearTarget} && -n ${yearSource} ]] \
    && path=${2/"/${yearTarget}"/"/${yearSource}"}
  echo ${path}
  return 0
}

guessYearFast()
{
  #guess the year from the path, pick the rightmost one
  #is string AUTOYEAR present, will be returned
  local IFS="/"
  declare -a pathArray=( ${1} )
  local field=""
  local year=""
  local autoYear=""
  for field in ${pathArray[@]}; do
    [[ ${field} =~ ^20[0-9][0-9]$ ]] && year="${field}"
    [[ ${field} == AUTOYEAR ]] && autoYear="${field}"
  done
  [[ -n ${autoYear} ]] && year="${autoYear}"
  echo ${year}
  return 0
}

run2year()
{
  #for a given run print the year.
  #the run-year table is ${PWGPP_runMap} (a string)
  #one line per year, format: year runMin runMax
  local run=$1
  [[ -z ${run} ]] && return 1
  local year=""
  local runMin=""
  local runMax=""
  while read year runMin runMax; do
    [[ -z ${year} || -z ${runMin} || -z ${runMax} ]] && continue
    [[ ${run} -ge ${runMin} && ${run} -le ${runMax} ]] && echo ${year} && break
  done < <(echo "${PWGPP_runMap}")
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
        echo "ALICE_ROOT=$ALICE_ROOT"
        echo "which aliroot: "$(which aliroot)
        echo "git describe:"
        echo "  "$(git -C $ALICE_ROOT/../src/ describe)
        echo 
        echo --------------------------------------

#
# ALIPHYSICS info
#
        echo --------------------------------------
        echo
        echo "ALIPHYSICSINFO"
        echo 
        echo "ALICE_PHYSICS=$ALICE_PHYSICS"
        echo "git describe:"
        echo "  "$(git -C $ALICE_PHYSICS/../src/ describe)
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

summarizeLogs()
{
  #validate and summarize the status of logs
  #input is a list of logs, or a glob:
  #example (summarizes logs in current and subdirs):
  #  summarizeLogs * */*
  #if no args given, process all files in PWD
  #exit code 1 if some logs are not validated

  #print a summary of logs
  local input
  local file=""
  declare -A files
  input=("${@}")
  [[ -z "${input[*]}" ]] && input=( "${PWD}"/* )
  
  #double inclusion protection+make full paths
  for file in "${input[@]}"; do
    [[ ! "${file}" =~ ^/ ]] && file="${PWD}/${file}"
    files["${file}"]="${file}"
  done

  local logFiles
  logFiles="\.*log$|^stdout$|^stderr$"

  #check logs
  local logStatus=0
  local errorSummary=""
  local validationStatus=""
  declare -A coreFiles
  for file in "${files[@]}"; do
    [[ ! -f ${file} ]] && continue
    #keep track of core files for later processing
    [[ "${file##*/}" =~ ^core$ ]] && coreFiles[${file}]="${file}" && continue
    [[ ! "${file##*/}" =~ ${logFiles} ]] && continue
    errorSummary=$(validateLog ${file})
    validationStatus=$?
    [[ validationStatus -ne 0 ]] && logStatus=1
    if [[ ${validationStatus} -eq 0 ]]; then 
      #in pretend mode randomly report an error in rec.log some cases
      echo "${file} OK"
    elif [[ ${validationStatus} -eq 1 ]]; then
      echo "${file} BAD ${errorSummary}"
    elif [[ ${validationStatus} -eq 2 ]]; then
      echo "${file} OK MWAH ${errorSummary}"
    fi
  done

  #report core files
  for x in "${coreFiles[@]}"; do
    echo ${x}
    chmod 644 ${x}
    #gdb --batch --quiet -ex "bt" -ex "quit" aliroot ${x} > stacktrace_${x//\//_}.log
    gdb --batch --quiet -ex "bt" -ex "quit" aliroot ${x} > stacktrace.log
    local nLines[2]
    #nLines=($(wc -l stacktrace_${x//\//_}.log))
    nLines=($(wc -l stacktrace.log))
    if [[ ${nLines[0]} -eq 0 ]]; then
      #rm stacktrace_${x//\//_}.log
      rm stacktrace.log
    else
      logStatus=1
      echo "${x%/*}/stacktrace.log"
    fi
  done

  return ${logStatus}
}

validateLog()
{
  #validate one log file
  #input is path to log file
  #output an error summary on stdout
  #exit code is 0 if validated, 1 otherwise
  log=${1}
  errorConditions=(
            'There was a crash'
            'floating'
            'error while loading shared libraries'
            'std::bad_alloc'
            's_err_syswatch_'
            'Thread [0-9]* (Thread'
            'AliFatal'
            '\.C.*error:.*\.h: No such file'
            'segmentation'
            'Segmentation fault'
            'Interpreter error recovered'
            ': command not found'
            ': comando non trovato'
            'core dumped'
  )

  warningConditions=(
            'This is serious'
  )

  local logStatus=0
  local errorSummary=""
  local warningSummary=""
  local errorCondition=""
  for errorCondition in "${errorConditions[@]}"; do
    local tmp=$(grep -m1 -e "${errorCondition}" ${log})
    local error=""
    [[ -n ${tmp} ]] && error=" : ${errorCondition}"
    errorSummary+=${error}
  done

  local warningCondition=""
  for warningCondition in "${warningConditions[@]}"; do
    local tmp=$(grep -m1 -e "${warningCondition}" ${log})
    local warning=""
    [[ -n ${tmp} ]] && warning=" : ${warningCondition}"
    warningSummary+=${warning}
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

mergeSysLogs()
{
  if [[ $# -lt 2 ]]; then
    echo 'merge syslogs to an output file'
    echo 'usage:'
    echo 'mergeSysLogs outputFile inputFile1 inputFile2 ...'
    return 0
  fi
  local outputFile
  local inputFiles
  local i
  local x
  local runNumber
  outputFile=${1}
  shift
  inputFiles="$@"
  i=0
  if ! ls -1 ${inputFiles} &>/dev/null; then echo "the files dont exist!: ${inputFiles}"; return 1; fi
  while read x; do 
    runNumber=$(guessRunNumber ${x})
    [[ -z ${runNumber} ]] && echo "run number cannot be guessed for ${x}" && continue
    awk -v run=${runNumber} -v i=${i} 'NR > 1 {print run" "$0} NR==1 && i==0 {print "run/I:"$0}' ${x}
    (( i++ ))
  done < <(ls -1 ${inputFiles}) > ${outputFile}
  return 0
}

stackTraceTree()
{
  if [[ $# -lt 1 ]]; then
    echo 'make stacktrace processing  in case of standard root crash log'
    echo 'input is a (list of) text files with the stack trace (either gdb aoutput'
    echo 'produced with e.g. gdb --batch --quiet -ex "bt" -ex "quit" aliroot core,'
    echo 'or the root crash log), output is a TTree formatted table.'
    echo 'example usage:'
    echo 'benchmark.sh stackTraceTree /foo/*/rec.log'
    echo 'benchmark.sh stackTraceTree $(cat file.list)'
    echo 'benchmark.sh stackTraceTree `cat file.list`'
    return 0
  fi
  #cat "${@}" | gawk '
  gawk '
       BEGIN { 
       print "frame/I:method/C:line/C:cpass/I:aliroot/I:file/C";
               RS="#[0-9]*";
               aliroot=0;
               read=1;
             } 
      /There was a crash/ {read=1;}
      /The lines below might hint at the cause of the crash/ {read=0;}
      read==1 { 
               if ($3 ~ /Ali*/) aliroot=1; else aliroot=0;
               gsub("#","",RT); 
               if ($NF!="" && RT!="" && $3!="") print RT" "$3" "$NF" "0" "aliroot" "FILENAME
             }
      ' "${@}" 2>/dev/null
}

plotStackTraceTree()
{
  #plot the stacktrace tree,
  #first arg    is the text file in the root tree format
  #second arg   is optional: a plot is written to file instead of screen
  #third arg    is optional: selection for plotting, default skip G_ stuff
  local tree=$1
  local plot=${2:-"crashes.png"}
  local selection=${3:-'!strstr(method,\"G__\")'}
  [[ ! -f ${tree} ]] && echo "plotStackTraceTree: no input file given" && return 1
  aliroot -b <<EOF
TTree* t=AliSysInfo::MakeTree("${tree}");
TCanvas* canvas = new TCanvas("QA crashes","QA crashes",1);
t->Draw("method","${selection}","");
canvas->SaveAs("${plot}");
.q
EOF
  return 0
}

#these functions encode strings to and from a space-less form
#use when spaces are not well handled (e.g. in arguments to 
#commands in makeflow files, etc.
encSpaces()(echo "${1// /@@@@}")
decSpaces()(echo "${1//@@@@/ }")

get_realpath() 
{
  if [[ $# -lt 1 ]]; then
    echo "print the full path of a file or directory, like \"readlink -f\" on linux"
    echo "Usage:"
    echo "  get_realpath <someFileOrDir>"
    return 0
  fi
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
  elif [[ -d "$1" ]]; then
    if cd "$1" &>/dev/null; then
      local tmppwd="$PWD"
      cd - &>/dev/null
      echo "$tmppwd"
      return 0
    else
      return 1
    fi
  else
    # file *cannot* exist
    return 1 # failure
  fi
  # reassemble realpath
  echo "$tmppwd"/"${1##*/}"
  return 0 # success
}

printLogStatistics()
{
  #this function processes the summary logs and prints some stats
  #relies on the summary log format produced by summarizeLogs()
  # - how many checked logs in total
  # - number of each type of problem
  # example usage:
  #   printLogStatistics */*.log
  [[ ! -f $1 ]] && return 1
  echo "log statistics from: ${1%/*}"
  #cat "${@}" | awk '
  awk '
  BEGIN {nOK=0; nCores=0; nStackTraces=0;}
  /\/core/ {nCores++}
  /\/stacktrace.log/ {nStackTraces++}
  /OK/ {nOK++; nLogs++;}
  /BAD/ {
    nLogs++
    err=""
    write=0
    for (i=3; i<=NF; i++)
    { 
      if ($i ~ /^\:$/) 
        write=1
      else
        write=0

      if (write==0)
      {
        if (err=="") err=$i
        else err=(err FS $i)
      }

      if (err != "" && (write==1 || i==NF))
      {
        sumBAD[err]++
        err=""
      }
    }
  } 
  END {
    print ("number of succesful jobs: " nOK" out of "nLogs )
    for (key in sumBAD)
    {
      print key": "sumBAD[key]
    }
    if (nCores>0) print "core files: "nCores", stack traces: "nStackTraces 
  }
  ' "${@}"
}

createUniquePID()
{
  #create a unique ID for jobs running in parallel
  #consists of the ip address of the default network interface, PID,
  #if an argument is given, append it (e.g. a production ID)
  #the fields are space separated with a tag for easy parsing
  #spaces in the productionID will be encoded using encSpaces()
  local productionID=""
  [[ -n "${1}" ]] && productionID=$(encSpaces "${1}")
  local defaultIP=$(/sbin/route | awk '$1=="default" {print $8}' | xargs /sbin/ifconfig | awk '/inet / {print $2}' | sed 's/.*\([0-9]?\.[0-9]?\.[0-9]?\.[0-9]?\)/$1/')
  local id="ip:${defaultIP} pid:${BASHPID}"
  [[ -n "${productionID}" ]] && id+=" prod:${productionID}"
  echo "${id}"
}

copyFileToLocal()
(
  #copies a single file to a local destination: the file may either come from
  #a local filesystem or from a remote location (whose protocol must be
  #supported)
  #copy is "robust" and it is repeated some times in case of failure before
  #giving up (1 is returned in that case)
  #origin: Dario Berzano, dario.berzano@cern.ch
  src="$1"
  dst="$2"
  ok=0
  [[ -z "${maxCopyTries}" ]] && maxCopyTries=10

  proto="${src%%://*}"

  echo "copy file to local dest started: $src -> $dst"

  for (( i=1 ; i<=maxCopyTries ; i++ )) ; do

    echo "...attempt $i of $maxCopyTries"
    rm -f "$dst"

    if [[ "$proto" == "$src" ]]; then
      cp "$src" "$dst"
    else
      case "$proto" in
        root)
          xrdcp -f "$src" "$dst"
        ;;
        http)
          curl -L "$src" -O "$dst"
        ;;
        *)
          echo "protocol not supported: $proto"
          return 2
        ;;
      esac
    fi

    if [ $? == 0 ] ; then
      ok=1
      break
    fi

  done

  if [[ "$ok" == 1 ]] ; then
    echo "copy file to local dest OK after $i attempt(s): $src -> $dst"
    return 0
  fi

  echo "copy file to local dest FAILED after $maxCopyTries attempt(s): $src -> $dst"
  return 1
)

paranoidCp()
(
  #recursively copy files and directories
  #if target is a directory - it must exist!
  #to avoid using find and the like as they kill
  #the performance on some cluster file systems
  #does not copy links to avoid problems
  sourceFiles=("${@}")
  destination="${sourceFiles[@]:(-1)}" #last element
   #check if we are not trying to copy to the same structure
  if [ $destination == $sourceFiles ] ; then  
      echo paranoidCp INFO skip
      echo paranoidCp INFO destination== $destination 
      echo paranoidCp INFO sourceFiles== $sourceFiles
      return 1; 
  fi

  unset sourceFiles[${#sourceFiles[@]}-1] #remove last element (dst)
  #[[ ! -f "${destination}" ]] 
  for src in "${sourceFiles[@]}"; do
    if [[ -f "${src}" && ! -h  "${src}" ]]; then
      paranoidCopyFile "${src}" "${destination}"
    elif [[ -d "${src}" && ! -h "${src}" ]]; then
      src="${src%/}"
      dst="${destination}/${src##*/}"
      mkdir -p "${dst}"
      paranoidCp "${src}"/* "${dst}"
    fi
  done
)

paranoidCopyFile()
(
  #copy a single file to a target in an existing dir
  #repeat a few times if copy fails
  #returns 1 on failure, 0 on success
  src=$(get_realpath "${1}")
  dst=$(get_realpath "${2}")
  [[ -d "${dst}" ]] && dst="${dst}/${src##*/}"
  #some sanity check
  [[ -z "${src}" ]] && echo "$1 does not exist" && return 1
  [[ -z "${dst}" ]] && echo "$2 does not exist" && return 1
  #check if we are not trying to copy to the same file
  [[ "${src}" == "${dst}" ]] && echo "$dst==$src, not copying" && return 0
  ok=0
  [[ -z "${maxCopyTries}" ]] && maxCopyTries=10

  echo "paranoid copy started: $src -> $dst"
  for (( i=1 ; i<=maxCopyTries ; i++ )) ; do

    echo "...attempt $i of $maxCopyTries"
    #rm -f "$dst"
    cp -a -n "$src" "$dst"

    cmp -s "$src" "$dst"
    if [ $? == 0 ] ; then
      ok=1
      break
    else
      rm "${dst}"
    fi

  done

  if [[ "$ok" == 1 ]] ; then
    echo "paranoid copy OK after $i attempt(s): $src -> $dst"
    return 0
  fi

  echo "paranoid copy FAILED after $maxCopyTries attempt(s): $src -> $dst"
  return 1
)

generPWD(){ 
  #
  # generate semirandom pwd using 2 keys 
  # Example usage:
  # generPWD  myserviceaccount10 key11
  key0=$1
  key1=$2
  heslo0=`md5sum <<< "$key0 $key1" | cut -c 1-16`
  heslo=`echo $heslo0 | cut -c 1-8| awk '{print toupper($0)}'`
  heslo=$heslo`echo $heslo0 | cut -c 8-15| awk '{print tolower($0)}'`%
  echo $heslo;
}

#this makes debugging easier:
#executes the command given as an argument in this environment
#use case:
#  bashdb utilities.sh summarizeLogs * */*
[[ $# != 0 ]] && eval "$@" || true
