#!/bin/bash
#
# need pubsub_size_tag program visible

usage="<simplefile input directory> <number of events> [<output directory> (optional, default is current)]"

if [[ "$1" == "-h" || "$1" == "-help" || "$1" == "--help" ]] ; then
  echo Usage: $0 $usage
  exit
fi

if [ `file $1 | awk '{print $2}'` != "directory" ]; then
  echo Error: $1 not a directory
  echo Usage: $0 $usage
  exit
fi

case $2 in 

  *[!0-9+-]*|?*[-+]*|""|-|+)
    echo Error: $2 not an integer
    echo Usage: $0 $usage
    exit
    ;;

  *)
    if [ ! $2 -ge 0 ] ; then
      echo Error: $2 not a positive integer
      echo Usage: $0 $usage
      exit
    fi
    ;;

esac


dir=`pwd`

if [ ! -z "$3" ] ; then
  mkdir $3
  dir=$3
fi

for (( i = 0 ; i < $2 ; i++ ));  do
  suffix=`EventStorageExtractor -storagetype simplefile -storagename $1 -eventnr $i -metadatadump | grep fDataType | grep fDataBlocks | awk '{print $4}' | sed s/\(// | tr A-Z a-z`
  EventStorageExtractor -storagetype simplefile -storagename $1 -eventnr $i -datadump > $dir/event_$i.$suffix
  ../bin/pubsub_size_tag -c $dir/event_$i.$suffix
done
