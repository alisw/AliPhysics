#!/bin/bash

# set job and simulation variables as :
# ./simrun.sh  --RUNNUMBER --NUMBEROFEVENTS --MONITOR(Bool_t)

function runcommand(){
    
    if [ "$2" = "run_mcgenFT2.C" ]; then
	echo "* $1 : $2($4,$5,$6,$7)"
	echo "* $1 : $2($4,$5,$6,$7)" >&2
	time aliroot -b -q -x $2\($4\,$5\,$6\,$7\) > $3 2>&1
    fi
    if [ "$2" = "CheckSmearedESD.C" ];then
	echo "* $1 : $2"
	echo "* $1 : $2" >&2
	time aliroot -b -q -x $2 > $3 2>&1
    fi
    if [ "$2" = "run_BAna.C" ];then
	echo "* $1 : $2"
	echo "* $1 : $2" >&2
	time aliroot -b -q -x $2\($4\) > $3 2>&1
    fi    
}

CONFIG_SEED=""
runnumber="$1"
events="$2"
monitor="$3"

CONFIG_SEED=$((ALIEN_PROC_ID%1000000000))

if [ "$CONFIG_SEED" -eq 0 ]; then
    CONFIG_SEED=$(((runnumber*100000+events)%1000000000))
    echo "* MC Seed is $CONFIG_SEED (based on run / counter : $runnumber / $monitor)"
else
    echo "* MC Seed is $CONFIG_SEED (based on AliEn job ID)"
fi

if [ "$CONFIG_SEED" -eq 0 ]; then
    echo "*!  WARNING! Seeding variable for MC is 0 !" >&2
fi

export ALIMDC_RAWDB1="./mdc1"
export ALIMDC_RAWDB2="./mdc2"
export ALIMDC_TAGDB="./mdc1/tag"
export ALIMDC_RUNDB="./mdc1/meta"

echo "SIMRUN:: Run $runnumber Event $events Monitor $monitor Seed $CONFIG_SEED"

runcommand "SIMULATION AND RECONSTRUCTION" "run_mcgenFT2.C" simRecAod.log $events $runnumber $monitor $CONFIG_SEED

runcommand "CHECK SMEARED ESD" "CheckSmearedESD.C" check.log

#runcommand "B Analysis" "run_BAna.C" Banalysis.log $runnumber
#mv syswatch.log completeWatch.log


exit 0
