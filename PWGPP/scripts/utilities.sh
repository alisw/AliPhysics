#!/usr/bin/env bash
#library of useful PWGPP related bash functions
#it REQUIRES BASH 4 !!!!
#blame: Mikolaj Krzewicki, mkrzewic@cern.ch

if [ ${BASH_VERSINFO} -lt 4 ]; then
  echo "bash version >= 4 needed, you have ${BASH_VERSION}, exiting..."
  exit 1
fi

# Load alilog4bash.sh from the same directory containing this script.
source "$(dirname "${BASH_SOURCE[0]}")"/alilog4bash.sh false

PWGPP_runMap="
2010 136833 139517 pbpb
2010 108350 136832 pp
2011 165747 170593 pbpb
2011 140441 165746 pp
2012 188230 188366 ppb
2012 171590 193766 pp
2013 195344 197388 ppb
2013 197469 197692 pp
2014 200008 208364 NONE
2015 244908 246994 pbpb
2015 224956 244628 pp
2016 249954 999999 pp
"

parseConfig()
{
  # parse command line arguments, they have to be in the form
  #  option=value
  #    --or--
  #  -option value
  #    --or--
  #  --option value
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
  local var=""
  local value=""
  for opt in "${args[@]}"; do
    [[ -n ${encodedSpaces} ]] && opt="$(decSpaces ${opt})"
    if [[ ! "${opt}" =~ .*=.* ]]; then
      if [[ "${opt}" =~ ^-.? && -z "$expectPosixOptionValue" ]]; then
        var="${opt#--}"
        var="${var#-}"
        expectPosixOptionValue=1
        continue
      elif [[ -n "$expectPosixOptionValue" ]]; then
        value="${opt}"
        unset expectPosixOptionValue
      else
        continue;  # non option string should be allowed - e.g parsing parameters for alihadd
        #echo "badly formatted option ${var}, should be: option=value (or -var value) stopping..."
        #return 1
      fi
    else
      var="${opt%%=*}"
      value="${opt#*=}"
    fi
    #echo "setting ${var}=${value}"
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
  #one line per year, format: year runMin runMax collisionSystem
  local run=$1
  [[ -z ${run} ]] && return 1
  local year=""
  local runMin=""
  local runMax=""
  local collisionSystem
  while read year runMin runMax collisionSystem; do
    [[ -z ${year} || -z ${runMin} || -z ${runMax} ]] && continue
    [[ ${run} -ge ${runMin} && ${run} -le ${runMax} ]] && echo ${year} && break
  done < <(echo "${PWGPP_runMap}")
  return 0
}

run2collisionSystem()
{
  #for a given run print the year.
  #the run-year table is ${PWGPP_runMap} (a string)
  #one line per year, format: year runMin runMax collisionSystem
  local run=$1
  [[ -z ${run} ]] && return 1
  local year=""
  local runMin=""
  local runMax=""
  while read year runMin runMax collisionSystem; do
    [[ -z ${year} || -z ${runMin} || -z ${runMax} ]] && continue
    [[ ${run} -ge ${runMin} && ${run} -le ${runMax} ]] && echo ${collisionSystem} && break
  done < <(echo "${PWGPP_runMap}")
  return 0
}

gitInfo(){
    #
    # print git information in alilog format - to enable parsing
    # in case rquested  diff file created in the $ALICE_ROOT and $ALICE_PHYSICS instalation directory
    #
    # USAGE:
    #     during code development to keep track of the software version  ( git describe ) as we can not use tags
    #
    makeDiff=$1
    if [[ -z ${makeDiff} ]] ; then
       echo "gitInfo <makeDiff>"; 
       makeDiff=2;
    fi
    alilog_info  "utilities.sh/gitInfo  START"
    alilog_info  "wdir="`pwd`
    alilog_info  "\$ALICE_ROOT="$ALICE_ROOT
    alilog_info  "\$ALICE_ROOT git describe="$(git -C $ALICE_ROOT/../src/ describe)
    alilog_info  "\$ALICE_PHYSICS="$ALICE_PHYSICS
    alilog_info  "\$ALICE_PHYSICS git describe="$(git -C $ALICE_PHYSICS/../src/ describe)
    alilog_info  "\$ROOTSYS="$ROOTSYS
    alilog_info  "\$ROOTSYS git describe="$(git -C $ROOTSYS/../src/ describe)
    if [ $makeDiff -eq 1 ] ; then    # dump diff file  to the install directory
	alilog_info "git  -C $ALICE_ROOT/../src/ diff >\$ALICE_ROOT/ALICE_ROOT.diff"
	git  -C $ALICE_ROOT/../src/ diff >$ALICE_ROOT/ALICE_ROOT.diff
	alilog_info "git  -C $ALICE_PHYSICS/../src/ diff >\$ALICE_PHYSICS/ALICE_PHYSICS.diff"
	git  -C $ALICE_PHYSICS/../src/ diff >$ALICE_PHYSICS/ALICE_PHYSICS.diff
    fi;
    if [ $makeDiff -eq 2 ] ; then      # copy software diff  if exist to the current directory
	alilog_info "cp -f $ALICE_ROOT/ALICE_ROOT.diff `pwd`"
	cp -f $ALICE_ROOT/ALICE_ROOT.diff .
	alilog_info "cp -f $ALICE_PHYSICS/ALICE_PHYSICS.diff `pwd`"
	cp -f $ALICE_PHYSICS/ALICE_PHYSICS.diff .
    fi;
    alilog_info  "utilities.sh/gitInfo  END"
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

parseListOfFiles()
{
  #generate a list of files, one per line out of arguments.
  #names starting with "@" are assumed to be file lists and will
  #be expanded
  for file in "$@"; do
    if [[ ${file:0:1} == @ ]]; then
      [[ -r "${file:1}" ]] && cat "${file:1}"
    else
      echo "$file"
    fi
  done
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
  local -a input
  local file=""
  declare -A files
  while IFS= read x; do
    input+=("$x")
  done < <(parseListOfFiles "$@")
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
    [[ ! -f "${file}" ]] && continue
    #keep track of core files for later processing
    [[ "${file##*/}" =~ ^core ]] && coreFiles[${file}]="${file}" && continue
    [[ ! "${file##*/}" =~ ${logFiles} ]] && continue
    errorSummary=$(validateLog "${file}")
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
    echo "core ${x}"
    chmod 644 "${x}"
    stacktraceLog=${x}.stacktrace.log
    #gdb --batch --quiet -ex "bt" -ex "quit" aliroot ${x} > stacktrace_${x//\//_}.log
    gdb --batch --quiet -ex "bt" -ex "quit" aliroot "${x}" > "$stacktraceLog"
    local nLines[2]
    #nLines=($(wc -l stacktrace_${x//\//_}.log))
    nLines=($(wc -l "$stacktraceLog"))
    if [[ ${nLines[0]} -eq 0 ]]; then
      #rm stacktrace_${x//\//_}.log
      rm "$stacktraceLog"
    else
      logStatus=1
      echo "stacktrace $stacktraceLog"
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
  log="${1}"
  [[ ! -f "$log" ]] && return 1
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
            'core dumped'
            'core file'
            'Exception catched'
            'line [0-9]*.*: No such file or directory'
            'core\.'
            'std::bad_alloc'
            'Segmentation violation'
            'Bus error'
            'floating point exception'
            'Killed'
            'busy flag cleared'
            'Cannot build the PAR archive'
            'glibc detected'
            'E-AliCDBGrid::PutEntry:'
            'F-AliCDBGrid::'
            'E-TAlienFile::ReadBuffer: The remote'
            'Compilation failed'
            'Low statistics: number of contributing tracks'
  )

  warningConditions=(
            'This is serious'
  )

  local logStatus=0
  local errorSummary=""
  local warningSummary=""
  local errorCondition=""
  for errorCondition in "${errorConditions[@]}"; do
    local tmp=$(grep -m1 -e "${errorCondition}" "${log}")
    local error=""
    [[ -n "${tmp}" ]] && error=" ; ${errorCondition}"
    errorSummary+=${error}
  done

  local warningCondition=""
  for warningCondition in "${warningConditions[@]}"; do
    local tmp=$(grep -m1 -e "${warningCondition}" "${log}")
    local warning=""
    [[ -n "${tmp}" ]] && warning=" : ${warningCondition}"
    warningSummary+="${warning}"
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
  if [[ $# -lt 1 ]]; then
    echo 'merge syslogs to an output file'
    echo 'usage:'
    echo 'mergeSysLogs outputFile inputFile1 inputFile2 ...'
    echo 'if file name prepended with "@" it is a file list'
    return 0
  fi
  local outputFile
  local inputFiles
  local i
  local x
  local runNumber
  local fsize
  outputFile=${1}
  shift
  inputFiles="$@"
  fileNumber=0
  parseListOfFiles ${inputFiles[@]} | while IFS= read logFile; do
    if [[ $fileNumber == 0 ]]; then
      tableHeader=$( head -n 1 "$logFile" )
      tableHeader+=":year/I:period/C:run/I:pass/C:line/D:itree/D:nstamps/D:treeName/C"
      echo $tableHeader
    fi

    guessRunData "$logFile"
    [[ -z "$year" ]] && year=0
    [[ -z "$period" ]] && period=0
    [[ -z $runNumber ]] && runNumber=0
    [[ -z "$pass" ]] && pass=0
    
    nLines=$(wc -l 2>/dev/null < "$logFile"); ((nLines--))

    lineNumber=-1
    while read line; do
      ((lineNumber++))
      [[ $lineNumber == 0 ]] && continue
      
      extraBranches="$year $period $runNumber $pass $lineNumber $fileNumber $nLines"
      echo "$line $extraBranches"

    done < "$logFile"
    
    (( fileNumber++ ))
  done > "$outputFile"
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
    echo 'benchmark.sh stackTraceTree @file.list somefile.log'
    return 0
  fi
  parseListOfFiles "$@" | while IFS= read x; do echo "filename: $x"; cat "$x" 2>/dev/null; done | gawk '
       BEGIN { 
       print "frame/I:method/C:line/C:cpass/I:aliroot/I:file/C";
               RS="#[0-9]*";
               aliroot=0;
               read=1;
             }
      /^filename:/ {filename=$2}
      /There was a crash/ {read=1;}
      /The lines below might hint at the cause of the crash/ {read=0;}
      read==1 { 
               if ($3 ~ /Ali*/) aliroot=1; else aliroot=0;
               gsub("#","",RT); 
               if ($NF!="" && RT!="" && $3!="") print RT" "$3" "$NF" "0" "aliroot" "filename
             }
      ' 2>/dev/null
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
  echo "log statistics from: ${1%/*}"
  parseListOfFiles "$@" | while IFS= read x; do cat "$x" 2>/dev/null; done | gawk '
  BEGIN {nOK=0; nCores=0; nStackTraces=0; nLogs=0;}
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
    print ("number of validated logs: " nOK" out of "nLogs )
    for (key in sumBAD)
    {
      print key": "sumBAD[key]
    }
    if (nCores>0 || nStackTraces>0) print "core files: "nCores", stack traces: "nStackTraces 
  }
  '
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

printExec() {
  alilog_info "[command] [PWD=$PWD]: $*" >&2
  "$@"
}

bigEcho() {
  # Outputs big text to stderr.
  which figlet &> /dev/null \
    && figlet -- "$@" >&2 \
    || ( printf "\n\n$*\n\n" | tr '[[:lower:]]' '[[:upper:]]' >&2 )
  true
}

listDir() (
  dir="$(cd "$1"; pwd)"
  echo
  alilog_info "[listDir] Content of ${dir}${2:+" ($2)"}"
  find "$dir" -ls
  echo
)

mkdirLocal() (
  # Creates the given directories, with full path, only if local. If not local,
  # print a message and return success.
  # Return 1 if at least one error happened, 0 on success.
  err=0
  while [[ $# -gt 0 ]]; do
    dir=$1
    shift
    if [[ "${dir%%://*}" != "$dir" ]]; then
      alilog_warning "[mkdirLocal] skipping creation of $dir: it is not local"
      continue
    fi
    mkdir -p "$dir"
    [[ -d "$dir" ]]
    rv=$?
    err=$((err + ($rv & 1)))
    [[ $rv == 0 ]] && alilog_success "[mkdirLocal] creation of dir $dir OK" \
                   || alilog_error   "[mkdirLocal] creation of dir $dir FAILED"
  done
  return $((err & 1))
)

statRemote() (
  # Check if file exists, whether it is local or remote. Returns 0 on success, 1 on failure.
  file=$1
  proto="${file%%://*}"
  [[ "$proto" == "$file" ]] && proto=local
  case "$proto" in
    local) [[ -f $file ]] ;;
    root)  path=${file:$((${#proto}+3))}
           host=${path%%/*}
           path=${path:$((${#host}))}
           while [[ ${path:0:2} == // ]]; do path=${path:1}; done
           xrd "$host" stat "$path" 2>&1 | grep -q Modtime: ;;
    *)     alilog_error "[statRemote] for path $file: protocol not supported: $proto"
           return 2 ;;
  esac
  rv=$?
  alilog_info "[statRemote] path $file (proto=${proto}$([[ $proto == local ]] && echo ", pwd=$PWD")) $([[ $rv == 0 ]] && echo "exists" || echo "does NOT exist")"
  return $rv
)

lsRemote() (
  # List remote files on stdout. Returns 0 on success, 1 on failure.
  dir=$1
  proto="${dir%%://*}"
  [[ "$proto" == "$dir" ]] && proto=local
  case "$proto" in
    local) find "$dir" -maxdepth 1 -mindepth 1 ;;
    root)  path=${dir:$((${#proto}+3))}
           host=${path%%/*}
           path=${path:$((${#host}))}
           while [[ ${path:0:2} == // ]]; do path=${path:1}; done
           xrd "$host" ls "$path" 2>&1 | grep -v "No such file or directory" | \
                                         grep -o "${path}.*" | \
                                         sed -e 's|^\(.*\)$|'$proto://$host/'\1|' ;;
    *)     alilog_error "[lsRemote] for directory $dir: protocol not supported: $proto" >&2
           return 2 ;;
  esac
  rv=$?
  [[ $rv != 0 ]] && alilog_error "[lsRemote] cannot list directory $dir" >&2
  return $rv
)

copyFileFromRemote() (
  # Copy a list of remote files ($1, $2...${n-1}) to a certain dest dir ($n).
  # The last parameter must then be a local directory. Files can be local too.
  # Dest dir is created if not existing.
  # If only one parameter is specified download this file to the current directory.
  # If a source file is in the form @list.txt, then the list of files (one per
  # line) in list.txt will be expanded and copied.
  # On success 0 is returned - 1 otherwise.
  # Example: copyFileFromRemote root://localhost//file1.txt @list.txt /tmp/localdir/foo

  [[ $# == 1 ]] && dstdir=. || dstdir=${!#}
  maxCopyTries=${maxCopyTries-10}
  sleepRetry=0
  err=0
  if which timeout &>/dev/null; then
    timeoutPrefix="timeout -s 9 ${remoteCpTimeout-600}"
  elif which gtimeout &>/dev/null; then
    timeoutPrefix="gtimeout -s 9 ${remoteCpTimeout-600}"
  else
    timeoutPrefix=""
  fi
  $timeoutPrefix true || timeoutPrefix=
  opname="[copyFileFromRemote]"
  mkdir -p "$dstdir"

  while [[ $# -gt 1 ]]; do
    [[ ${1:0:1} == @ ]] && inputcmd="cat ${1:1}" || inputcmd="echo $1"
    while read -u 4 src; do
      thiserr=1
      proto="${src%%://*}"
      [[ "$proto" == "$src" ]] && proto=local
      dst="$dstdir/$(basename "$src")"
      alilog_info "$opname (proto=$proto) started: $src -> $dst"
      for ((i=1; i<=maxCopyTries; i++)); do
        [[ $i -gt 1 ]] && alilog_warning "$opname $src -> $dst failed $i out of $maxCopyTries time(s), retrying in $sleepRetry s"
        sleep $sleepRetry
        sleepRetry=$((sleepRetry+1))
        case "$proto" in
          local) fullsrc=$(get_realpath "$src")
                 fulldst=$(get_realpath "$dst")
                 printExec $timeoutPrefix cp "$src" "$dst"
                 if [[ $? != 0  && "$fullsrc" != "$fulldst" ]]; then
                   rm -f "$dst"
                   false
                 fi ;;
          root)  printExec $timeoutPrefix xrdcp -f "$src" "$dst" ;;
          http*) printExec $timeoutPrefix curl -LsSfo "$dst" "$src" ;;
          *)     alilog_error "protocol not supported: $proto"
                 return 2 ;;
        esac
        if [[ $? == 0 ]]; then
          thiserr=0
          break
        fi
      done  # cp attempts
      [[ $thiserr == 0 ]] && alilog_success "$opname (proto=$proto) OK after $i attempt(s): $src -> $dst" \
                          || alilog_error   "$opname (proto=$proto) FAILED after $maxCopyTries attempt(s): $src -> $dst"
      err=$((err+thiserr))
    done 4< <($inputcmd)
    shift
  done

  return $((err & 1))
)

copyFileToRemote() (
  # Copy a list of files ($1, $2...${n-1}) to a certain dest dir ($n). The last
  # parameter must then be a directory. The destination might be local or remote
  # if a protocol is specified (e.g. root://).
  # Dest dir is created if not existing (xrdcp does that already).
  # If a source file is in the form @list.txt, then the list of files (one per
  # line) in list.txt will be expanded and copied.
  # Shell glob is supported.
  # On success 0 is returned - 1 otherwise.
  # Example: copyFileToRemote file1.txt files* @list.txt test1 root://host/dir

  dstdir=${!#}
  maxCopyTries=${maxCopyTries-10}
  proto="${dstdir%%://*}"
  sleepRetry=0
  err=0
  if which timeout &>/dev/null; then
    timeoutPrefix="timeout -s 9 ${remoteCpTimeout-600}"
  elif which gtimeout &>/dev/null; then
    timeoutPrefix="gtimeout -s 9 ${remoteCpTimeout-600}"
  else
    timeoutPrefix=""
  fi
  $timeoutPrefix true || timeoutPrefix=
  [[ "$proto" == "$dstdir" ]] && proto=local
  opname="[copyFileToRemote] (proto=$proto)"

  while [[ $# -gt 1 ]]; do
    [[ ${1:0:1} == @ ]] && inputcmd="cat ${1:1}" || inputcmd="echo $1"
    while read -u 4 src; do
      thiserr=1
      dst="$dstdir/$(basename "$src")"
      for ((i=1; i<=maxCopyTries; i++)); do
        [[ -d "$src" ]] && echo "$opname $src -> $dst skipping, is a directory" && thiserr=0 && break
        # TODO use shopt -s nullglob
        [[ ! -r "$src" ]] && echo "$opname $src -> $dst skipping, cannot access source" && thiserr=0 && break
        [[ $i -gt 1 ]] && alilog_warning "$opname $src -> $dst failed $i out of $maxCopyTries time(s), retrying in $sleepRetry s"
        sleep $sleepRetry
        sleepRetry=$((sleepRetry+1))
        case "$proto" in
          local) mkdir -p "$(dirname "$dst")"
                 fullsrc=$(get_realpath "$src")
                 fulldst=$(get_realpath "$dst")
                 printExec $timeoutPrefix cp "$src" "$dst"
                 if [[ $? != 0 && "$fullsrc" != "$fulldst" ]]; then
                   rm -f "$dst"
                   false
                 fi ;;
          root)  printExec $timeoutPrefix xrdcp -f "$src" "$dst" ;;
          *)     alilog_error "protocol not supported: $proto"
                 return 2 ;;
        esac
        if [[ $? == 0 ]]; then
          thiserr=0
          break
        fi
      done  # cp attempts
      [[ $thiserr == 0 ]] && alilog_success "$opname OK after $i attempt(s): $src -> $dst" \
                          || alilog_error   "$opname FAILED after $maxCopyTries attempt(s): $src -> $dst"
      err=$((err+thiserr))
    done 4< <($inputcmd)
    shift
  done

  return $((err & 1))
)

function xCopy() {
  # Recursive and parallel copy for files. Usage:
  #   xCopy -w 10 -d proto://host//destpath file1 file2 @list dir...
  # Where:
  #   -d is the destination directory (can be local, or remote)
  #   -w is the number of parallel workers (defaults to 15)
  #   -f does not preserve the source dir structure and copies all flat to destpath
  #   -c skips copy if local source does not exist
  #   -C skips copy if local destination already exists
  # Note that if one of the sources is a local directory, it will be inspected recursively. This
  # does not apply to remote directories instead.

  opname="[xCopy]"
  workers=0
  flat=0
  checklocalsrc=0
  checklocaldest=0
  while [[ "$1" != -- ]]; do
    case "$1" in
      --destdir|-d) dstdir="$2"; shift 2 ;;
      --flat|-f) flat=1; shift ;;
      --check-local-src|-c) checklocalsrc=1; shift ;;
      --check-local-dest|-C) checklocaldest=1; shift ;;
      --workers|-w) workers=$(($2)); shift 2 ;;
      *) break ;;
    esac
  done
  [[ "$dstdir" == '' ]] && { alilog_error "$opname Missing --destdir" ; exit 1 ; }
  [[ $workers == 0 ]] && workers=15  # default
  alilog_info "$opname Copying files to $dstdir using $workers workers"
  [[ "${dstdir%%://*}" == "$dstdir" ]] && copyFunc=copyFileFromRemote || copyFunc=copyFileToRemote

  # Create a directory of symlinks: readlink will tell us what is each source
  # file, with respect to the *current* directory. Symlinks are used because
  # operations on them are atomic.
  t=$(mktemp -d /tmp/xcp.XXXXX)
  count=0
  while [[ $# -gt 0 ]]; do
    [[ ${1:0:1} == @ ]] && inputcmd="cat ${1:1}" || inputcmd="echo $1"
    while read src; do
      [[ "${src%%://*}" == "$src" && -d "$src" ]] && inputcmd="find $src -type f" \
                                                  || inputcmd="echo $src"
      while read -u 4 src2; do
        [[ $checklocalsrc == 1 && "${src2%%://*}" == "$src2" && ! -f "${src2}" ]] \
          && { alilog_warning "Skipping local nonexisting source $src2" ; continue ; }
        ln -nfs $src2 $t/$count
        count=$((count+1))
      done 4< <($inputcmd)
    done < <($inputcmd)
    shift
  done

  # Start workers.
  for ((i=0; i<workers; i++)); do
    ( while [[ 1 ]]; do
        placeholder=$(find $t -type l -print -quit 2> /dev/null)
        [[ "$placeholder" == '' ]] && break
        src=$(readlink $placeholder)
        rm $placeholder 2> /dev/null || continue
        dstdir_flat=$dstdir/$([[ $flat == 0 ]] && dirname $src)
        [[ $copyFunc == copyFileFromRemote && $checklocaldest == 1 \
                                           && -f $dstdir_flat/$(basename $src) ]] \
          && { alilog_warning "Skipping $src: local destination exists" ; continue ; }
        $copyFunc $src $dstdir_flat
      done
    ) &
  done

  # Do not leave rubbish behind if dying. Kills the whole process group (-$$).
  trap "rm -rf $t; kill -9 -$$" SIGHUP SIGINT SIGTERM

  wait
  rm -rf $t
}

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
      rm -f "${dst}"
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

reformatXMLCollection()
{
  #parse the xml collection on stdinput
  #output a collection in format: file md5 ctime size
  local nfields=""
  local turl=""
  local md5=""
  local ctime=""
  local size=""
  local x=""
  while read -a fields; do
    nfields=${#fields[*]}
    turl=""
    md5=""
    ctime=""
    size=""
    for ((x=1;x<=${nfields};x++)); do
      field=${fields[${x}]}
      if [[ "${field}" == "md5="* ]]; then
        eval ${field}
      fi
      if [[ "${field}" == "turl="* ]]; then
        eval ${field}
      fi
      if [[ "${field}" == "ctime="* ]]; then
        eval "${field} ${fields[((x+1))]}"
      fi
      if [[ "${field}" == "size="* ]]; then
        eval ${field}" "${fields[((x+1))]}
      fi
    done
    ctime=$( date -d "${ctime}" +%s 2>/dev/null)
    [[ -z $md5 ]] && md5="."
    [[ -n "$turl" ]] && echo "${turl//"alien://"/} ${md5} ${ctime} ${size}"
  done
}

reportDoneFile()
{
  #print a report on a file if exists
  #args: tag filename directory
  local tag=$1
  local file=$2
  local dir=$3
  [[ -r $file ]] && echo "$tag ${dir}/${file}"
}


mergeAliSysInfo(){
    #
    # marian.ivanov@cern.ch
    # I think this method should be in $ALICE_ROOT/STEER/AliSysInfo.sh script  
    # 
    # mergeAliSysInfo trees
    # in addition to the original branches file information is appended to the tree
    # New info appended
    # Parameters:
    #   1.) inputList  as ascii file
    #   2.) outputFile 
    #       in case outputFile is root file - tree is written as root file therwise plane txt used
    # example usage 
    # ( source $ALICE_PHYSICS/../src/PWGPP/scripts/utilities.sh;  mergeAliSysInfo syswatchHis.list syswatchHis.root; )
    # ( source $ALICE_PHYSICS/../src/PWGPP/scripts/utilities.sh;  mergeAliSysInfo syswatchMap.list syswatchMap.root; )
    inputList=$1
    outputFile=$2
    alilog_info "mergeAliSysInfo  inputList=$1 outputFile=$2 BEGIN"
    if [[ -z ${inputList} ]] ; then
       echo "mergeAliSysInfo [inputList] [outputFile]"; 
       return 0;
    fi
    if [[ -z ${outputFile} ]] ; then
       echo "mergeAliSysInfo [inputList] [outputFile]"; 
       return 0;
    fi
    #
    counter=0
    desc=`head  -n 1 $inputList  | xargs head -n 1`  
    echo $desc:line/D:itree/D:nstamps/D:treeName/C  >$outputFile.tree
    fsize=0
    for afile in `cat  $inputList`; do 
	fsize=$(wc -l < $afile)
	((counter++))
	#echo $counter $fsize $afile 
	lineCounter=0
	cat $afile | grep -v "hname" | while read -r line; 	
	  do 
	  ((lineCounter++))
	  echo "$line" $lineCounter $counter $fsize $afile; 
	done 
    done >> $outputFile.tree
    [[ $outputFile =~ .root$ ]] &&  echo "AliSysInfo::MakeTree(\"$outputFile.tree\",\"$outputFile\")" | aliroot -b 
    alilog_info "mergeAliSysInfo  inputList=$1 outputFile=$2 END"
    return 1;
}





#this makes debugging easier:
#executes the command given as an argument in this environment
#use case:
#  bashdb utilities.sh summarizeLogs * */*
[[ $# != 0 ]] && eval "$@" || true


# vi:syntax=zsh
