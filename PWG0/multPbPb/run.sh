#!/bin/bash

# defaults
isMC=kFALSE
run=no
correct=no
nev=-1
offset=0
debug=kFALSE
runmode=1
dataset=/PWG3/zampolli/run000104867_90_92_pass4
ropt="-l"
option="SAVE"
workers=26
analysismode=9; #SPD + field on

give_help() {

cat <<ENDOFGUIDE
This scripts runs the mupliplicity analysis.

Available options:
 Mode control, at least one of the following options should be used
  -r <mode>                    Run the task
                               Modes [default=$runmode]:
                                  0 local
                                  1 caf    
  -c                           Run the correction
 Proof settings
  -w nworkers                  Set the number of worker nodes
  -n <nev>                     Number of events to be analized 
 Misc
  -d <dataset>                 Dataset or data collection (according to run mode) [default=$dataset]
  -o <option>                  Misc option [default=$option]
                               Available options: 
                                - SAVE:  move results to a different output folder
                                - ITSsa: Use ITSsa tracks
                                - TPC:   Use TPC only tracks
  -t <option>                  Command line option for root [defaul=$ropt]
  -m                           Use this to run on Monte Carlo
  -g                           Debug mode
  -h                           This help
ENDOFGUIDE

}

while getopts "r:cgmd:o:w:n:" opt; do
  case $opt in
    r)
      run=yes
      runmode=$OPTARG
      ;;      
    c)
      correct=yes
      ;;      
    g)
      debug=kTRUE;
      ;;
    m)
      isMC=kTRUE
      ;;
    d)
      dataset=$OPTARG
     ;;
    o)
      option=$OPTARG
      ;;
    t)
      ropt=$OPTARG
      ;;
    w) 
      workers=$OPTARG
      ;;
    n) 
      nev=$OPTARG
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

if [ "$run" = "$correct" ]
    then 
    echo "One and only one option between -r and -c must be selected"
    give_help
    exit 1
fi

if [ "$run" = "yes" ]
    then
    root $ropt run.C\(\"$dataset\",$nev,$offset,$debug,$runmode,$isMC,\"$option\",$workers\)
fi
