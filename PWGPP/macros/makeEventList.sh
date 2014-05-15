#!/bin/bash
main()
{
  #create the event list for one run: the run number if guessed from the path
  #the input file can either be a root file with the trees or a zipfile with 
  #file FilterEvents_Trees.root inside
  [[ $# -eq 0 ]] && echo "Usage: $0 file" && exit
  file=$1
  outputFile=$2

  runNumber=$(guessRunNumber $file)
  period=$(guessPeriod $file)
  ptMinHighPt="8.0";
  ptMinV0s="3.0";
  if [[ ${period%_*} =~ (LHC10h|LHC11h|LHC12h) ]]; then
    ptMinHighPt="14.0";
    ptMinV0s="6.0";
  fi
  
  [[ -z $outputFile ]] && outputFile=${runNumber}.list

  [[ ! -f $file ]] && echo "cannot access file $file" && exit

  [[ "${file##*\.}" == *zip ]] && file+="#FilterEvents_Trees.root"

  this=$0
  [[ ${OSTYPE} =~ inux ]] && this=$(readlink -f $0)
  scriptPath=${this%/*}

  echo outputFile=$outputFile
  echo runNumber=$runNumber
  echo period=$period
  echo PWD=$PWD
  echo scriptPath=$scriptPath
  echo aliroot -b -q "${scriptPath}/makeEventList.C(\"${file}\",${ptMinHighPt},${ptMinV0s})"

  aliroot -b -q "${scriptPath}/makeEventList.C(\"${file}\",${ptMinHighPt},${ptMinV0s})" 2>/dev/null\
  | awk -v period=${period} '/^offlineTrigger/ {triggerType=$2;} $8 !~ "*" && $0 ~ "^*.*/\\w*/\\w*/" { n=split($10,a,"/"); rawfile="/"a[4]"/"a[5]"/"a[6]"/"period"/"a[8]"/raw/"a[n-1]".root"; print rawfile" "$8" "triggerType; }' \
  | sort -V | uniq > $outputFile
}

guessRunNumber()
{
  (
  #guess run number from the path, pick the rightmost one
  IFS="/"
  declare -a path=( $1 )
  dirDepth=${#path[*]}
  for ((x=${dirDepth}-1;x>=0;x--)); do
    field=${path[${x}]}
    [[ ${field} =~ ^000[0-9]*[0-9]$ ]] && runNumber=${field#000} && break
  done
  echo $runNumber
  )
}

guessPeriod()
{
  (
  #guess the year from the path, pick the rightmost one
  IFS="/"
  declare -a path=( $1 )
  dirDepth=${#path[*]}
  for ((x=${dirDepth}-1;x>=0;x--)); do
    field=${path[${x}]}
    [[ ${field} =~ ^LHC[0-9][0-9][a-z]_*.*$ ]] && period=${field%%_*} && break
  done
  echo $period
  )
}

main "$@"
