#!/bin/bash

# Example usage:
#
#   source $ALICE_ROOT/libexec/alihadd.sh
#   # Filter input list before merging
#   makeFilteredList input.list highPt filter.log
#   # Internal function
#   testFileList input.list highPt | tee  filter.log

loadUtilities() {
  #  load utilities if not yet done
  if [[ -z $(type -t alilog_info) ]]; then
    source $ALICE_ROOT/libexec/alilog4bash.sh
    alilog_info "Load ALICE_ROOT/libexec/alilog4bash.sh"
  fi
  if [[ -z $(type -t guessPeriod) ]]; then
    source $ALICE_ROOT/libexec/utilities.sh
    alilog_info "Load $ALICE_ROOT/libexec/utilities.sh"
  fi
}

testFileList() {
  # Internal function: filter input list
  # Provide good/bad file list out of inputList(param[$1])
  # Tests: read/write/delete/eventually copytree for selected keys(param[$2])
  inputList=$1
  keyList=$2
  timeOut=$3
  [[ $timeOut > 0 ]] && prefix="timeout $timeOut" || prefix=
  alilog_info "inputList=$inputList"
  alilog_info "keyList=$keyList"
  for a in $(cat $inputList); do
    alilog_info "Filter.Start $a $keyList $timeOut"
    ( export ROOT_HIST=0
      $prefix aliroot -n -b -l <<EOF
AliXRDPROOFtoolkit::TestFile("$a", "$keyList");
EOF
    )
    alilog_info "Filter.End $a  $?"
  done
}

processFilterLog() {
  # Internal function: parse log file produced by testFileList function
  inputList=$1
  logFile=$2
  alilog_info "processFilterLog $logFile"
  # get stat
  ninput=$(cat $logFile | grep -c "Filter.Start ")
  nall=$(cat $logFile | grep  -c "I-AliXRDPROOFtoolkit::TestFile.TestAll:.*TestStatusAll_*")
  ngood=$(cat $logFile | grep  -c -e "I-AliXRDPROOFtoolkit::TestFile.TestAll:.*TestStatusAll_0")
  alilog_info "processFilterLog $ninput:$nall:$ngood"
  cat $logFile | grep "TestStatusAll_0" | awk '{print $2}' > $inputList.Good
  cat $logFile | grep "TestStatusAll_" | grep -v "TestStatusAll_0" | awk '{print $2}' > $inputList.Bad
}

makeFilteredList() {
  if [ $# -lt 3 ]; then
    printf "\n"
    alilog_error "makeFilteredList.InvalidParameters $@"
    printf "\n"
    printf "Usage: makeFilteredList inputList keyList logfile <timeOut>\n"
    return
  fi
  inputList=$1
  keyList=$2
  logFile=$3
  timeOut=$4
  testFileList $inputList $keyList $timeOut | tee $logFile
  processFilterLog $inputList $logFile
}

alihaddsh() {
  # the same as alihadd but additional flags
  #
  # Question  - hot to parse flags
  # - checkKeys pattern
  # - prefetch prefix
  # source $ALICE_ROOT/libexec/alihadd.sh
  # alihaddsh -k -s 1000000000   -checkKeys highPt -prefetch 1 -v 1 -timeOut 300 test.root @input.list
  # alihaddsh -k -s 1000000000   -checkKeys highPt -prefetch 1 -v 1 -timeOut 300 test.root -inputList input.list
  alilog_info "alihaddsh $*"
  config="$@"
  parseConfig $config
  # check @inputList in parameter string if not specified as   -inputList
  [[ -z "$inputList" ]] && inputList=$(echo $config|grep " @"| sed -e 's/.*@\(.*\)/\1/'| gawk '{print $1}')
  if [[ -z "$inputList" ]]; then
    alilog_error "Input list not specified"
    exit 1
  else
    alilog_info "Input list $inputList"
  fi

  if [[ -n checkKeys ]]; then
    alilog_info "alihaddsh checkKeys $checkKeys"
    config=$(echo $config |sed s_-checkKeys.*-_-_)
    if [ $timeOut  > 0 ]; then
      alilog_info "alihaddsh timeOut $timeOut"
      config=$(echo $config |sed s_-timeOut\ [0-9]*\ __)
    fi
    alilog_info "makeFilteredList $inputList $checkKeys filter.log $timeOut"
    makeFilteredList $inputList $checkKeys filter.log $timeOut
    config=$(echo $config | sed s/${inputList}/${inputList}.Good/)
  fi
  if [[ -n prefetch ]]; then
    alilog_info "alihaddsh prefetch $prefetch"
    alilog_info "NOT YET IMPLEMENTED"
    #config=$(echo $config |sed s_-prefetch.*-_-_)
  fi
  alilog_info "Starting: alihadd $config"
  alihadd $config
  if [[ $? == 0 ]]; then
    alilog_success "alihadd.sh $config"
    return 0
  else
    alilog_error "alihadd.sh $config"
    return 1
  fi
}

loadUtilities
eval "$@"
