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
runTriggerStudy=no
customSuffix=""
ntrackletsTrigger=50
rejectBGV0Trigger=kFALSE
useTrackCentralityCut=0
trackMin=0
trackMax=100

give_help() {

cat <<ENDOFGUIDE
This scripts runs the mupliplicity analysis and the trigger study task

Available options:
 Mode control, at least one of the following options should be used
  -r <mode>                    Run the task
                               Modes [default=$runmode]:
                                  0 local
                                  1 caf    
                                  2 grid    
  -c                           Run the correction
  -s                           Run the trigger study task (by default it runs the multiplicity analysis)
 Proof settings
  -w nworkers                  Set the number of worker nodes
  -n <nev>                     Number of events to be analized 
 Misc
  -d <dataset>                 Dataset or data collection (according to run mode) [default=$dataset]
                                - local mode: a single ESD file, an xml collection of files on 
                                  grid or a text file with a ESD per line
                                - caf mode: a dataset
                                - grid mode: a directory on alien
 Options specific to the multiplicity analysis
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
  -y <min,max>                 Select centrality based on "good tracks" rather than on centrality
                               estimator [off by default]
  -0 <min,max>                 Select centrality based on v0 multiplicity range rather than on centrality
                               estimator [off by default]
  -o <option>                  Misc option [default=$option]
                               Available options: 
                                - SAVE:     Move results to a different output folder*
                                - ITSsa:    Use ITSsa tracks
                                - TPC:      Use TPC only tracks
                                - NOMCKINE: Skip MC kinematics (runs way faster)
                                * == can be used in trigger studies task
  -t <option>                  Command line option for root [defaul=$ropt]
  -m                           Use this to run on Monte Carlo
  -x  <suffix>                 Set a custom suffix in the histo manager        
  -g                           Debug mode
  -h                           This help

 Options specific to the trigger study task
  -k <ntracklets>              Max number of tracklets to fill eta and pt 
                               distributions [default=$ntrackletsTrigger]
  -v                           Reject BG with the V0
ENDOFGUIDE

}

while getopts "x:sr:cgmd:o:w:n:e:b:t:k:vy:0:" opt; do
  case $opt in
    r)
      run=yes
      runmode=$OPTARG
      ;;      
    y)
      useTrackCentralityCut=1
      trackMin=${OPTARG%%,*}
      trackMax=${OPTARG##*,}
      ;;      
    0)
      useTrackCentralityCut=2
      trackMin=${OPTARG%%,*}
      trackMax=${OPTARG##*,}
      ;;      
    x)
      customSuffix=$OPTARG
      ;;      
    k)
      ntrackletsTrigger=$OPTARG
      ;;      
    s)
      runTriggerStudy=yes
      ;;      
    v)
      rejectBGV0Trigger=kTRUE
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
    if [ "$runTriggerStudy" = "yes" ]
	then
	root $ropt runTriggerStudy.C\(\"$dataset\",$nev,$offset,$debug,$runmode,$isMC,$ntrackletsTrigger,$rejectBGV0Trigger,\"$option\",$workers\)
    else
	root $ropt run.C\(\"$dataset\",$nev,$offset,$debug,$runmode,$isMC,$centrBin,\"$centrEstimator\",$useTrackCentralityCut,$trackMin,$trackMax,\"$option\",\"$customSuffix\",$workers\)
    fi
fi

if [ "$correct" = "yes" ]
    then
    echo "To be implemented"
fi
