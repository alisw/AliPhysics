#!/bin/bash

run=137161
pass=pass1_5plus
mc=0
mode="full"
nev=1234566789
workers=26
ROPT=""
listfile=""
offset=0
debug=kTRUE
option="SAVE" #FIXME:set option
suffix=""

give_help() {

cat <<ENDOFGUIDE
This scripts runs the the physics selection and centrality on a specified run/dataset

Available options:
 Mode control, at least one of the following options should be used
  -r <mode>                    Run the task
                               Modes (compulsory):
                                  0 local
                                  1 caf
                                  2 grid    
  -d <run or dataset>          Run number(s) (grid) or dataset (caf) or file name (local)
                               Local filename can be an xml collection of files on alies, 
                               a single esd, or a text file with an ESD filename per line
  -m                           Set if runnning on MC
  -t <rootopt>                 Options passed to root
  -l <list.txt>                Process sequentially all runs/dataset listed in the file 
                               (one entry per line). If you use this option you don't 
                               need the -d. in the case you are running on CAF, they 
                               must have the same path
 Grid only options
  -g <gridmode>                Plugin Mode [default=$mode]
  -p <recopass>                Reconstruction pass [default=$pass]       
 CAF only options
  -p <path>                    Data set path
  -n <nev>                     Number of events
  -w <workers>                 Number of workers [default=$workers]
ENDOFGUIDE

}

while getopts "r:hd:mg:p:n:w:t:l:" opt; do
  case $opt in
    r)
      runMode=$OPTARG
      ;;
    d)
      run=$OPTARG
      ;;
    l)
      listfile=$OPTARG
      ;;
    t)
      ROPT=$OPTARG
      ;;
    n)
      nev=$OPTARG
      ;;
    w)
      workers=$OPTARG
      ;;
    m)
      mc=kTRUE
      ;;
    g) 
      mode=$OPTARG
      ;;
    p)
      pass=$OPTARG
      ;;
    h)
      give_help
      exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      give_help
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      give_help
      exit 1
      ;;
  esac
done


runlist=$run
if [ "$listfile" != "" ]
    then
    runlist=""
    while read line
    do
	runlist="$runlist $line"
    done < $listfile
    
fi

echo "Run list: $runlist"


if [ "$runMode" = "2" ]
then
    echo root $ROPT run.C\(\"$run\",\"$pass\",$nev,$offset,$debug,$runMode,$mc,$option,$suffix,$workers,\"$mode\"\)
    root $ROPT run.C\(\"$run\",\"$pass\",$nev,$offset,$debug,$runMode,$mc,\"$option\",\"$suffix\",$workers,\"$mode\"\)
else
    for run in $runlist 
    do
	echo root $ROPT run.C\(\"$run\",\"$pass\",$nev,$offset,$debug,$runMode,$mc,$option,$suffix,$workers,\"$mode\"\)
	root $ROPT run.C\(\"$run\",\"$pass\",$nev,$offset,$debug,$runMode,$mc,\"$option\",\"$suffix\",$workers,\"$mode\"\)
    done
fi
