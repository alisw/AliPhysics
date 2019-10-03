#!/bin/bash

Usage() {
    echo 'Usage: mrgAlgOut.sh <source_dir1> <source_dir2> ...'
    echo 'Merged output will be stored in source directory'
    exit 1
}

control=mpControlRes.root
stat=mpStatOut.root
mpdat=mpData.root

# determine own directory
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$SDIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

if [ $# -lt 1 ] ; then Usage ;fi

for dir in "$@"; 
do
#
 dir=`echo $dir | sed "s/\/$//"`
 echo doing $dir
 ${SDIR}/mrgfl.sh -f ${dir}/$control ${dir}/*/$control
 ${SDIR}/mrgfl.sh -f ${dir}/$stat ${dir}/*/$stat
 ${SDIR}/mrgfl.sh -f ${dir}/$mpdat ${dir}/*/$mpdat

#
done

