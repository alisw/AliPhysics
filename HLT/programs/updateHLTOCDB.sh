#!/usr/bin/env bash
ocdbSource="/cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/"
ocdbTarget="OCDB/2015/"
cdbEntries="*/*/*"
firstRun=""

main ()
{
  if [[ -z "$@" ]]; then
    echo "update the OCDB dir with the most recent versions of selected CDB objects"
    echo "extend the validity of those objects to infinity"
    echo "normally the target CDB will already be filled with the default entries"
    echo "and only the entries without defaults (some GRP/*/* entries) will be missing"
    echo "* first prepare a copy of the default objects using e.g.:"
    echo '* aliroot -b -q $ALICE_SRC'"'"'/HLT/programs/downloadCDB.C(999999,"local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB","local://OCDB/2015","*/*/*")'"'"
    echo
    echo "example:"
    echo "  ${0##*/} ocdbSource=/cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/ ocdbTarget=/OCDB/2015/ cdbEntries=*/*/*"
    return 1
  fi

  if ! parseConfig "$@"; then
    ${0}
    return 1
  fi

  for sourceEntryPath in ${ocdbSource}/${cdbEntries}
  do 
    cdbEntry=${sourceEntryPath#$ocdbSource/}
    targetEntryPath=${ocdbTarget}/${cdbEntry}

    #latestFileTarget=$(/bin/ls -1 ${targetEntryPath} | sort --field-separator=_ --key=2Vr,3Vr | grep -v Run0_999 | head -1)
    #defaulFileTarget=$(/bin/ls -1 ${targetEntryPath} | sort --field-separator=_ --key=2Vr,3Vr | grep    Run0_999 | head -1)
    #runRangeLatestTarget=( $(getOCDBrunRange ${latestFileTarget##*/}) )
    #runRangeDefaultTarget=( $(getOCDBrunRange ${defaulFileTarget##*/}) )
    
    latestFile=$(/bin/ls -1 ${sourceEntryPath} | sort --field-separator=_ --key=2Vr,3Vr | grep -v Run0_999 | head -1)
    #defaulFile=$(/bin/ls -1 ${sourceEntryPath} | sort --field-separator=_ --key=2Vr,3Vr | grep    Run0_999 | head -1)
    runRangeLatest=( $(getOCDBrunRange ${latestFile##*/}) )
    #runRangeDefault=( $(getOCDBrunRange ${defaulFile##*/}) )
    
    #if this is empty it means we only have the default object (or none)
    #in which case we do nothing, because it should already be copied
    [[ -z ${runRangeLatest[1]} ]] && continue

    if [[ -z ${firstRun} ]]; then
      runRange=( "${runRangeLatest[@]}" )
    else
      runRange=( ${firstRun} "999999999" "v0" "s0" )
    fi

    #TODO: this is stupid
    if [[ ! -d ${targetEntryPath} ]]; then
      firstRun=${runRange[0]}
      copyOCDBentry ${cdbEntry} local://${ocdbSource} local://${ocdbTarget} ${firstRun} "999999999"
    fi
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
