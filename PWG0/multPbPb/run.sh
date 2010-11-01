#!/bin/bash

# defaults
isMC=kFALSE
run=no
correct=no
nev=-1
offset=0
debug=kFALSE
runmode=1
dataset=/alice/sim/LHC10f8a_130844
ropt="-l"
option="SAVE"
workers=26
analysismode=9; #SPD + field on
centrBin=0
centrEstimator="V0M"

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
  -b <bin>                     Set centrality bin [default=$centrBin]
  -e <estimator>               Set centrality estimator [default=$centrEstimator]
                               Available choiches:
                                - V0M = V0 multiplicity
                                - FMD = FMD raw multiplicity
                                - TRK = N. of tracks
                                - TKL = N. of tracklets
                                - CL0 = N. of clusters in layer 0
                                - V0MvsFMD = correlation between V0 and FMD
                                - TKLvsV0 = correlation between tracklets and V0
                                - ZEMvsZDC = correlation between ZEM and ZDC     
  -o <option>                  Misc option [default=$option]
                               Available options: 
                                - SAVE:     Move results to a different output folder
                                - ITSsa:    Use ITSsa tracks
                                - TPC:      Use TPC only tracks
                                - NOMCKINE: Skip MC kinematics (runs way faster)
  -t <option>                  Command line option for root [defaul=$ropt]
  -m                           Use this to run on Monte Carlo
  -g                           Debug mode
  -h                           This help
ENDOFGUIDE

}

while getopts "r:cgmd:o:w:n:e:b:" opt; do
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
    e)
      centrEstimator=$OPTARG
     ;;
    b)
      centrBin=$OPTARG
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
    root $ropt run.C\(\"$dataset\",$nev,$offset,$debug,$runmode,$isMC,$centrBin,\"$centrEstimator\",\"$option\",$workers\)
fi

if [ "$correct" = "yes" ]
    then
    echo "To be implemented"
fi
