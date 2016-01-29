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
    [[ "${file##*/}" =~ ^core$ ]] && coreFiles[${file}]="${file}" && continue
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
    echo "${x}"
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
            ': comando non trovato'
            'core dumped'
            'core file'
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
    [[ -n "${tmp}" ]] && error=" : ${errorCondition}"
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
  if [[ $# -lt 2 ]]; then
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
  outputFile=${1}
  shift
  inputFiles="$@"
  i=0
  parseListOfFiles "$inputFiles" | while IFS= read x; do
    runNumber=$(guessRunNumber "${x}")
    [[ -z ${runNumber} ]] && echo "run number cannot be guessed for ${x}" && continue
    gawk -v run=${runNumber} -v i=${i} 'NR > 1 {print run" "$0} NR==1 && i==0 {print "run/I:"$0}' "${x}"
    (( i++ ))
  done > "${outputFile}"
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
  echo "[command] [PWD=$PWD]: $*" >&2
  "$@"
}

listDir() (
  dir="$(cd "$1"; pwd)"
  echo ; echo "[listDir] Content of ${dir}${2:+" ($2)"}"
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
      echo "[mkdirLocal] skipping creation of $dir: it is not local"
      continue
    fi
    mkdir -p "$dir"
    [[ -d "$dir" ]]
    rv=$?
    err=$((err + ($rv & 1)))
    echo "[mkdirLocal] creation of dir $dir $([[ $rv == 0 ]] && echo "OK" || echo "FAILED")"
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
    *)     echo "[statRemote] for file $file: protocol not supported: $proto"
           return 2 ;;
  esac
  rv=$?
  echo "[statRemote] file $file (proto=${proto}$([[ $proto == local ]] && echo ", pwd=$PWD"))" \
       "$([[ $rv == 0 ]] && echo "exists" || echo "does NOT exist")"
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
           xrd "$host" ls "$path" 2>&1 | grep -v "No such file or directory" | grep -o "${path}.*" ;;
    *)     echo "[lsRemote] for directory $dir: protocol not supported: $proto" >&2
           return 2 ;;
  esac
  rv=$?
  [[ $rv != 0 ]] && echo "[lsRemote] cannot list directory $dir" >&2
  return $rv
)

copyFileFromRemote() (
  # Copy a list of remote files ($1, $2...${n-1}) to a certain dest dir ($n).
  # The last parameter must then be a local directory. Files can be local too.
  # Dest dir is created if not existing.
  # If a source file is in the form @list.txt, then the list of files (one per
  # line) in list.txt will be expanded and copied.
  # On success 0 is returned - 1 otherwise.
  # Example: copyFileFromRemote root://localhost//file1.txt @list.txt /tmp/localdir/foo

  dstdir=${!#}
  maxCopyTries=${maxCopyTries-10}
  err=0
  opname="[copyFileFromRemote]"
  mkdir -p "$dstdir"

  while [[ $# -gt 1 ]]; do
    [[ ${1:0:1} == @ ]] && inputcmd="cat ${1:1}" || inputcmd="echo $1"
    while read -u 3 src; do
      thiserr=1
      proto="${src%%://*}"
      [[ "$proto" == "$src" ]] && proto=local
      dst="$dstdir/$(basename "$src")"
      echo "$opname (proto=$proto) started: $src -> $dst"
      for ((i=1; i<=maxCopyTries; i++)); do
        echo "$opname $src -> $dst attempt $i of $maxCopyTries"
        case "$proto" in
          local) printExec cp "$src" "$dst"
                 if [[ $? != 0 ]]; then
                   rm -f "$dst"
                   false
                 fi ;;
          root)  printExec xrdcp -f "$src" "$dst" ;;
          http*) printExec curl -LsSfo "$dst" "$src" ;;
          *)     echo "protocol not supported: $proto"
                 return 2 ;;
        esac
        if [[ $? == 0 ]]; then
          thiserr=0
          break
        fi
      done  # cp attempts
      [[ $thiserr == 0 ]] && echo "$opname (proto=$proto) OK after $i attempt(s): $src -> $dst" \
                          || echo "$opname (proto=$proto) FAILED after $maxCopyTries attempt(s): $src -> $dst"
      err=$((err+thiserr))
    done 3< <($inputcmd)
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
  err=0
  [[ "$proto" == "$dstdir" ]] && proto=local
  opname="[copyFileToRemote] (proto=$proto)"

  while [[ $# -gt 1 ]]; do
    [[ ${1:0:1} == @ ]] && inputcmd="cat ${1:1}" || inputcmd="echo $1"
    while read -u 3 src; do
      thiserr=1
      dst="$dstdir/$(basename "$src")"
      for ((i=1; i<=maxCopyTries; i++)); do
        [[ -d "$src" ]] && echo "$opname $src -> $dst skipping, is a directory" && thiserr=0 && break
        # TODO use shopt -s nullglob
        [[ ! -r "$src" ]] && echo "$opname $src -> $dst skipping, cannot access source" && thiserr=0 && break
        echo "$opname $src -> $dst attempt $i of $maxCopyTries"
        case "$proto" in
          local) mkdir -p "$(dirname "$dst")"
                 printExec cp "$src" "$dst"
                 if [[ $? != 0 ]]; then
                   rm -f "$dst"
                   false
                 fi ;;
          root)  printExec xrdcp -f "$src" "$dst" ;;
          *)     echo "protocol not supported: $proto"
                 return 2 ;;
        esac
        if [[ $? == 0 ]]; then
          thiserr=0
          break
        fi
      done  # cp attempts
      [[ $thiserr == 0 ]] && echo "$opname OK after $i attempt(s): $src -> $dst" \
                          || echo "$opname FAILED after $maxCopyTries attempt(s): $src -> $dst"
      err=$((err+thiserr))
    done 3< <($inputcmd)
    shift
  done

  return $((err & 1))
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

#this makes debugging easier:
#executes the command given as an argument in this environment
#use case:
#  bashdb utilities.sh summarizeLogs * */*
[[ $# != 0 ]] && eval "$@" || true


# vi:syntax=zsh
