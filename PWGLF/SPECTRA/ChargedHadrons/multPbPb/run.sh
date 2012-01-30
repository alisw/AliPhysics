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
option="DCA,SAVE"
workers=26
analysismode=9; #SPD + field on
centrBin=-1
centrEstimator="V0M"
runTriggerStudy=no
customSuffix=""
ntrackletsTrigger=50
rejectBGV0Trigger=kFALSE
useTrackCentralityCut=0
trackMin=0
trackMax=100
dataDir=""
mcDir=""
vzMin=-10
vzMax=10
etaMin=-0.5
etaMax=0.5
npart=381.188
weakFactor=-1
useSingleBin=kTRUE
#OUTPATH=output.BAK2010
OUTPATH=output

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
  -c <data,mc>                 Run the correction data and MC are names of the folders. 
                               ./output/ is added automatically in front of the folder names
  -s                           Run the trigger study task (by default it runs the multiplicity analysis)
 Proof settings
  -w nworkers                  Set the number of worker nodes (0 == 1 worker per node)
  -n <nev>                     Number of events to be analized 
 Misc
  -d <dataset>                 Dataset or data collection (according to run mode) [default=$dataset]
                                - local mode: a single ESD file, an xml collection of files on 
                                  grid or a text file with a ESD per line
                                - caf mode: a dataset
                                - grid mode: a directory on alien
  -h                           This help
 Options specific to the multiplicity analysis
  -l                           Run over all centrality bins
  -o <option>                  Misc option [default=$option]
                               Available options: 
                                - SAVE:     Move results to a different output folder*
                                - DCA:      Use DCA cut with global tracks
                                - ITSsa:    Use ITSsa tracks
                                - TPC:      Use TPC only tracks
                                - NOMCKINE: Skip MC kinematics (runs way faster)
                                * == can be used in trigger studies task
  -t <option>                  Command line option for root [defaul=$ropt]
  -m                           Use this to run on Monte Carlo
  -x  <suffix>                 Set a custom suffix in the histo manager        
  -g                           Debug mode
  == The following options are only valid if running on a single bin ==
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
  -0 <min,max>                 Select centrality based on v0 amplitude range rather than on centrality
                               estimator [off by default]
  -2 <min,max>                 Select centrality based on SPD outer layer clusters  rather than on centrality
                               estimator [off by default]

 Options specific to the trigger study task
  -k <ntracklets>              Max number of tracklets to fill eta and pt 
                               distributions [default=$ntrackletsTrigger]
  -v                           Reject BG with the V0
 Options specific for the corrections
  -z <zmin,zmax>               Change vertex Z range [default = $vzMin,$vzMax]
  -a <etamin,etamax>           Change eta range [default = $etaMin,$etaMax]
  -p <npart>                   Number of participants, used only for dNdeta/npart [default=$npart]
  -k <weakFrac>                Scale ration secondaries from strangeness/all rec by this factor [default=$weakFactor]
  -b <bin>                     Set centrality bin to be corrected. Only valid if you processed multiple 
                               bins at one (it changes the suffix of the multPbPbtracks.root file). It it's -1,
                               a file without suffix is searched for. This options applyies both to the data and to the 
                               MC file. [default=$centrBin]
ENDOFGUIDE

}

while getopts "x:sr:c:gmd:o:w:n:e:b:t:k:vy:0:2:hz:a:lp:" opt; do
  case $opt in
    r)
      run=yes
      runmode=$OPTARG
      ;;      
    l)
      useSingleBin=kFALSE
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
    2)
      useTrackCentralityCut=3
      trackMin=${OPTARG%%,*}
      trackMax=${OPTARG##*,}
      ;;      
    x)
      customSuffix=$OPTARG
      ;;      
    p)
      npart=$OPTARG
      ;;      
    k)
      ntrackletsTrigger=$OPTARG
      weakFactor=$OPTARG
      ;;      
    s)
      runTriggerStudy=yes
      ;;      
    v)
      rejectBGV0Trigger=kTRUE
      ;;      
    c)
      correct=yes
      dataDir="./$OUTPATH/${OPTARG%%,*}"
      mcDir="./$OUTPATH/${OPTARG##*,}"
      ;;      
    z)
      vzMin=${OPTARG%%,*}
      vzMax=${OPTARG##*,}
      ;;      
    a)
      etaMin=${OPTARG%%,*}
      etaMax=${OPTARG##*,}
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
	root $ropt run.C\(\"$dataset\",$nev,$offset,$debug,$runmode,$isMC,$centrBin,\"$centrEstimator\",$useTrackCentralityCut,$trackMin,$trackMax,\"$option\",\"$customSuffix\",$workers,$useSingleBin\)
    fi
fi

if [ "$correct" = "yes" ]
    then
    root $ropt correct.C+\(\"$dataDir\",\"$mcDir\",$vzMin,$vzMax,$etaMin,$etaMax,$npart,$weakFactor,$centrBin\);
fi
