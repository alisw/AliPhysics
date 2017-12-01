#!/usr/bin/env bash
ocdbSource="HCDB_test1"
ocdbTarget="HCDB_test1_fixed"
cdbEntries="*/*/*"
firstRun="0"
lastRun="999999999"
sourceRun="241544"

main()
{
  if [[ -z "$@" ]]; then
    echo "extend/modify the validity of OCDB objects"
    echo "* first prepare a copy of the default objects for someRunNumbr using e.g.:"
    echo '* aliroot -b -q $ALICE_SRC'"'"'/HLT/programs/downloadCDB.C(someRunNumber,"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB","local://OCDB/2015","*/*/*")'"'"
    echo
    echo "example:"
    echo "  ${0##*/} ocdbSource=/cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/ ocdbTarget=OCDB/2015/ cdbEntries=*/*/* sourceRun=<someRunNumber> firstRun=<begin target validity range> lastRun=<end target validity range>"
    return 1
  fi

  if ! parseConfig "$@"; then
    ${0}
    return 1
  fi

  for sourceEntryPath in ${ocdbSource}/${cdbEntries}
  do 
    cdbEntry=${sourceEntryPath#$ocdbSource}
    targetEntryPath=${ocdbTarget}/${cdbEntry}

    copyOCDBentry ${cdbEntry} local://${ocdbSource} local://${ocdbTarget} ${firstRun} ${lastRun} ${sourceRun}
  done
}

getOCDBrunRange()
{
  #output runStart, runEnd, v, s
  local data=${1#Run}
  data=${data%.root}
  echo ${data//_/ }
}

copyOCDBentry()
{
  #copy the CDB entry to a new storage
  #also modify the run range
  if [ $# -eq 0 ]; then
    echo "copyOCDBentry cdbEntry ocdbSource ocdbTarget firstRun lastRun sourceRun"
    echo "copyOCDBentry /GRP/GRP/Data local:///cvmfs/.../ local://./OCDB 0 999999999 241412"
    return
  fi

  local cdbEntry="${1}"
  local ocdbSource="${2}"
  local ocdbTarget="${3}"
  local firstRun="${4}"
  local lastRun="${5-999999999}"
  local sourceRun="${6-0}"
  aliroot -b << EOF
man=AliCDBManager::Instance();
man->SetDefaultStorage("${ocdbSource}");

AliCDBEntry* entry = man->Get("${cdbEntry}",${sourceRun});

AliCDBId entryID = entry->GetId();
entryID.SetRunRange(${firstRun},${lastRun});
entry->SetId(entryID);

AliCDBStorage * targetStorage = AliCDBManager::Instance()->GetStorage("${ocdbTarget}");
printf("storing ${cdbEntry} with range: %i-%i in ${ocdbTarget}\n",${firstRun},${lastRun});
targetStorage->Put(entry);
EOF
}

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

main "$@"
