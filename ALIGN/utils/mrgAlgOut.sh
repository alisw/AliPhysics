#!/bin/bash

Usage() {
    echo 'Usage: mrgAlgOut.sh <source_dir1> <source_dir2> ...'
    echo 'Merged output will be stored in source directory'
    exit 1
}

if [ $# -lt 1 ] ; then Usage ;fi

for dir in "$@"; 
do
#
 dir=`echo $dir | sed "s/\/$//"`
 echo doing $dir
 hadd -f ${dir}/controlRes.root ${dir}/*/controlRes.root
 hadd -f ${dir}/statOut.root ${dir}/*/statOut.root
 hadd -f ${dir}/mpData.root ${dir}/*/mpData.root
#
done

