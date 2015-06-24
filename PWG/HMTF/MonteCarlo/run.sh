#!/bin/bash

# defaults

ALIEN_USERNAME=`alien_whoami`
ALIEN_USERNAME=`echo $ALIEN_USERNAME  | tr -d '[[:space:]]'` #trim whitespaces

nev=1000
dir=dummy
tune=kPyTunePerugia2011
#tune=kPyTuneMonash2013
energy=7000
COMMAND=cp
sendOnGrid=no
runSim="yes"
runTask="yes"
taskOnly="no"
runOnGrid="no"
outsuffix=""
etamax=0.5
HEPMCFILE="pythia.hepmc"
#obsolete: threshold="50,100,150,200,250,300,350,400,450,500,550,600"
#obsolete: scaling="52819,45406,33136,16998,9373,5701,3590,1819,641,132,17,1"

#Perugia2011 scaling factors:
threshold="50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800"
scaling="1463100,1162700,843300,527300,380300,272900,181500,104800,51200,20900,7100,2000,500,100,10,1"
#threshold="50,100,150,200,250,300,350,400,450,500,550,600,650,700"
#scaling="14631,11627,8433,5273,3803,2729,1815,1048,512,209,71,20,5,1"
#new up to 600:
#threshold="50,100,150,200,250,300,350,400,450,500,550,600"
#scaling="575,457,331,207,150,107,71,41,20,8,3,1"
#OBSOLETE
#scaling="572,454,329,206,149,107,70,41,20,8,3,1"

#no scaling:
#threshold="10,20"
#scaling="1,1"

LIBSOPTIONS="0";

rm -fr dummy

# non settable variables
TARGETDIR="alien:/alice/cern.ch/user/${ALIEN_USERNAME:0:1}/${ALIEN_USERNAME}/HMTF/Test/"
CURRENTDIR=`pwd`


###
give_help() {

    echo "Available Options: "
    echo " -d <dir>             Target dir (DEF: $dir)";
    echo " -t <tune>            Set PYTHIA Tune or phojet flag (DEF: $tune)";
    echo "                       - kPyTuneMonash2013 (Pythia 8.205 with a Monash 2013 tune)";
    echo "                       - kPyTuneCDFA (Tune A)";
    echo "                       - kPyTuneAtlasCSC (Atlas Tune)";
    echo "                       - kPyTuneCMS6D6T  (CDF tune used by CMS - aka CMS tune)";
    echo "                       - kPyTunePerugia0 (Perugia 0 tune)";
    echo "                       - kPyTunePerugia2011 (Perugia 2011)"
    echo "                       - kPhojet (run phojet rather than pythia)";
    echo "                       - kHepMC  (run the HepMC gen reader rather than pythia)";
    echo " -n <nev>             number of events to generate (DEF: $nev)";
    echo " -e <energy>          Set energy in GeV (DEF: $energy)";
    echo " -l                   link, rather than copying, macros to the target dir)";
    echo " -r <task/sim/both>   run only task, only sim or both (DEF: both)";
    echo " -g                   use this flag to send the job on the grid";
    echo " -b                   use this flag when the job is running on the grid "
    echo "                       (it should only be used as an argument in the JDL) ";
    echo " -s <suff>            Add suffix to .log files (DEF: \"$outsuffix\")";
    echo " -x <etamax>          dNdeta is computed in the range |etamax| < 0.5 (DEF: \"$etamax\")";
    echo " -i <file>            HepMC input file. Needed if tune == kHepMC (DEF: \"$HEPMCFILE\")"
    echo ""
    echo " This script can downscale multiplicity bins so that the mult distribution looks flatter. "
    echo " The scaling is controlled by the 2 parameters belos"
    echo " -m <multiplicity>    array of multiplicity thresholds (DEF: $threshold)";   
    echo " -f <factor>          array of scaling factors (DEF: $scaling)";
}

 
while getopts ":d:t:n:m:f:e:r:s:x:hlgbi:" opt; do
  case $opt in
    d)
      dir=$OPTARG
      ;;
    i)
      HEPMCFILE=$OPTARG
      ;;
    t)
      tune=$OPTARG
#      if [ $tune = "phojet" ]
#	  then
#	  echo "Using phojet"
#	  tune="Tune_"
#      fi
      ;;
    r)
	  arg=$OPTARG
	  if [ $arg = "task" ]
	      then
	      runTask="yes"
	      runSim="no"
	      taskOnly="yes"
	  elif [ $arg = "sim" ]
	      then
	      runTask="no"
	      runSim="yes"
	  else
	      runTask="yes"
	      runSim="yes"
	  fi
     ;;
    n)
      nev=$OPTARG
      ;;
    m)
      threshold="$OPTARG"
      ;;
    f)
      scaling="$OPTARG"
      ;;
    x)
      etamax=$OPTARG
      ;;
    g)
      sendOnGrid="yes"
      ;;
    b)
      runOnGrid="yes"
      ;;
    s)
      outsuffix=$OPTARG
      ;;
    e)
      energy=$OPTARG
      ;;
    l)
      COMMAND="ln -s"
      ;;
    h) give_help
       exit
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done



if [ $runOnGrid != "yes" -a $sendOnGrid != "yes" -a $taskOnly != "yes" ]
    then
    if [ -e $dir ]
	then
	echo "Dir already exists. please remove it and try again"
        exit
	echo "Dir $dir already exists, remove (yes/no)?";
	read ans
	if [ "$ans" = "yes" ]
	    then
	    rm -rf $dir
	else
	    exit;
	fi
    fi    
    mkdir -p output/$dir
    ln -s output/$dir
    cd $dir
    $COMMAND ${CURRENTDIR}/LoadLibs.C .
    $COMMAND ${CURRENTDIR}/rungen.C .
    $COMMAND ${CURRENTDIR}/runTask.C .
    $COMMAND ${CURRENTDIR}/AliAnalysisTaskHMTFMC.cxx .
    $COMMAND ${CURRENTDIR}/AliAnalysisTaskHMTFMC.h .
fi

if [ $taskOnly = "yes" ] 
    then
    #previous if was skipped, but we need a fresh copy of task and runTask macro
    cd output/$dir
    $COMMAND ${CURRENTDIR}/runTask.C .
    $COMMAND ${CURRENTDIR}/AliAnalysisTaskHMTFMC.cxx .
    $COMMAND ${CURRENTDIR}/AliAnalysisTaskHMTFMC.h .

fi

echo "GRID $sendOnGrid"

if [ $sendOnGrid = "yes" ] 
    then

    pwd
    #remove old files
    alien_rm ${TARGETDIR/alien:/}/LoadLibs.C 
    alien_rm ${TARGETDIR/alien:/}/rungen.C 
    alien_rm /alice/cern.ch/user/${ALIEN_USERNAME:0:1}/${ALIEN_USERNAME}/bin/runFastSimulation.sh
    alien_rm ${TARGETDIR/alien:/}/simdndeta.jdl
    alien_rm ${TARGETDIR/alien:/}/validation.sh
    alien_mkdir ${TARGETDIR/alien:/}/$dir
    
    #upload new files
    # you can set an alien_CLOSE_SE if needed
    alien_cp  LoadLibs.C "$TARGETDIR"
    alien_cp  rungen.C "$TARGETDIR"
    alien_cp  validation.sh "$TARGETDIR"
    alien_cp  run.sh alien:/alice/cern.ch/user/${ALIEN_USERNAME:0:1}/${ALIEN_USERNAME}/bin/runFastSimulation.sh

    echo "USER: /${ALIEN_USERNAME:0:1}/${ALIEN_USERNAME}"

    sed s/MYARGS/"-d $dir -t $tune -n $nev -e $energy -m \"$threshold\" -f \"$scaling\" -b"/ $CURRENTDIR/templateJOB.jdl > tmp.jdl
    DIRJOB=${TARGETDIR//\//\\\/}
    DIRJOB=${DIRJOB/alien:/}
    sed s/TARGETDIR/$DIRJOB/ tmp.jdl > tmp2.jdl
    sed s/USERPATH/"${ALIEN_USERNAME:0:1}\/${ALIEN_USERNAME}"/ tmp2.jdl > tmp.jdl
    sed s/USER/"${ALIEN_USERNAME}"/ tmp.jdl > tmp2.jdl
    sed s/OUTDIR/$dir/ tmp2.jdl > simdndeta.jdl
    alien_cp simdndeta.jdl "$TARGETDIR"
else
    echo "WD: `pwd`"
    if [ "$runSim" = "yes" ]
    then
        if [ "$tune" = "kPhojet" ]
        then
           LIBSOPTIONS="-1"
        fi
        if [ "$tune" = "kPyTuneMonash2013" ]
        then
           LIBSOPTIONS="14"
        fi
	echo "GEN"
	echo root -b -q ./LoadLibs.C\($LIBSOPTIONS\) rungen.C+\($tune,$energy,$nev,\"$threshold\"\,\"$scaling\",\"$HEPMCFILE\"\)
        root -b -q ./LoadLibs.C\($LIBSOPTIONS\) rungen.C+\($tune,$energy,$nev,\"$threshold\"\,\"$scaling\",\"$HEPMCFILE\"\)  2>&1 | tee gen$outsuffix.log
        ls >> gen$outsuffix.log
    fi
    if [ "$runTask" = "yes" ]
	then
	echo "TASK"
	echo root -b -q ./LoadLibs.C\($LIBSOPTIONS\) runTask.C\($etamax\)
	root -b -q ./LoadLibs.C\($LIBSOPTIONS\) runTask.C\($etamax\) 2>&1 | tee task$outsuffix.log
    fi
fi

