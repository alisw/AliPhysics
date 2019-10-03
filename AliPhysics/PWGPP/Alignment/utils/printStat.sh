#!/bin/bash

Usage() {
    echo 'Usage: printStat.sh <mpStatOut1.root> <mpStatOut2.root> ...'
    exit 1
}


# determine own directory
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$SDIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

if [ $# -lt 1 ] ; then Usage ;fi

for stf in "$@"; 
do
#
# echo doing $stf
 aliroot -b -q -l ${SDIR}/printStat.C\(\"${stf}\"\) 
#
done

