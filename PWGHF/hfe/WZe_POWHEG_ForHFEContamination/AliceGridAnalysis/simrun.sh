#!/bin/bash

# set job and simulation variables as :
# ./simrun.sh  --run <x> --event <y> --process <kPythia6/kPhojet/kPythia6ATLAS_Flat/kPythia6D6T> --field <kNoField/k5kG> --energy <900/2360/10000>

function runcommand(){
    echo -e "\n"
    echo -e "\n" >&2

    echo "* $1 : $2"
    echo "* $1 : $2" >&2
    
    time aliroot -b -q -x $2 >>$3 2>&1
    exitcode=$?
    
    expectedCode=${5-0}
    
    if [ "$exitcode" -ne "$expectedCode" ]; then
        echo "*! $2 failed with exitcode $exitcode, expecting $expectedCode"
        echo "*! $2 failed with exitcode $exitcode, expecting $expectedCode" >&2
        echo "$2 failed with exitcode $exitcode, expecting $expectedCode" > validation_error.message
        exit ${4-$exitcode}
    else
        echo "* $2 finished with the expected exit code ($expectedCode), moving on"
        echo "* $2 finished with the expected exit code ($expectedCode), moving on" >&2
    fi
}

# Define the pt hard bin arrays
pthardbin_loweredges=( 0 5 11 21 36 57 84 117 152 191 234 )
pthardbin_higheredges=( 5 11 21 36 57 84 117 152 191 234 -1)

CONFIG_SEED=""
CONFIG_RUN_TYPE=""
CONFIG_FIELD=""
CONFIG_ENERGY=""
CONFIG_PHYSICSLIST=""
CONFIG_BMIN=""
CONFIG_BMAX=""
CONFIG_PTHARDBIN=""
CONFIG_PTHARDMIN=""
CONFIG_PTHARDMAX=""
CONFIG_QUENCHING=""
DC_RUN=""
DC_EVENT=""

RUNMODE=""

while [ ! -z "$1" ]; do
    option="$1"
    shift
    
    if [ "$option" = "--run" ]; then
            DC_RUN="$1"
            shift
    elif [ "$option" = "--event" ]; then
            DC_EVENT="$1"
            shift
    elif [ "$option" = "--process" ]; then
            CONFIG_RUN_TYPE="$1"
            shift
    elif [ "$option" = "--field" ]; then
            CONFIG_FIELD="$1"
            shift
    elif [ "$option" = "--energy" ]; then
            CONFIG_ENERGY="$1"
            shift
    elif [ "$option" = "--physicslist" ]; then
            CONFIG_PHYSICSLIST="$1"
            shift
    elif [ "$option" = "--bmin" ]; then
            CONFIG_BMIN="$1"
            shift
    elif [ "$option" = "--bmax" ]; then
            CONFIG_BMAX="$1"
            shift
    elif [ "$option" = "--pthardbin" ]; then
            CONFIG_PTHARDBIN="$1"
            shift  
    elif [ "$option" = "--quench" ]; then
            CONFIG_QUENCHING="$1"
            shift 
    elif [ "$option" = "--sdd" ]; then
            RUNMODE="SDD"
    fi
done

CONFIG_SEED=$((ALIEN_PROC_ID%1000000000))

if [ "$CONFIG_SEED" -eq 0 ]; then
    CONFIG_SEED=$(((DC_RUN*100000+DC_EVENT)%1000000000))
    echo "* MC Seed is $CONFIG_SEED (based on run / counter : $DC_RUN / $DC_EVENT)"
else
    echo "* MC Seed is $CONFIG_SEED (based on AliEn job ID)"
fi

if [ "$CONFIG_SEED" -eq 0 ]; then
    echo "*!  WARNING! Seeding variable for MC is 0 !" >&2
fi

echo "* b min is $CONFIG_BMIN"
echo "* b max is $CONFIG_BMAX"
echo "* pt hard bin is $CONFIG_PTHARDBIN"

if [ ! -z "$CONFIG_PTHARDBIN" ]; then
    # Define environmental vars for pt binning
    CONFIG_PTHARDMIN=${pthardbin_loweredges[$CONFIG_PTHARDBIN]}
    CONFIG_PTHARDMAX=${pthardbin_higheredges[$CONFIG_PTHARDBIN]}

    echo "* pt hard from $CONFIG_PTHARDMIN to $CONFIG_PTHARDMAX"
fi

mkdir input
mv galice.root ./input/galice.root
mv Kinematics.root ./input/Kinematics.root
ls input

export CONFIG_SEED CONFIG_RUN_TYPE CONFIG_FIELD CONFIG_ENERGY CONFIG_PHYSICSLIST CONFIG_BMIN CONFIG_BMAX CONFIG_PTHARDBIN CONFIG_PTHARDMIN CONFIG_PTHARDMAX DC_RUN DC_EVENT

export ALIMDC_RAWDB1="./mdc1"
export ALIMDC_RAWDB2="./mdc2"
export ALIMDC_TAGDB="./mdc1/tag"
export ALIMDC_RUNDB="./mdc1/meta"

if [ -f "$G4INSTALL/bin/geant4.sh" ]; then
    echo "* Sourcing G4 environment from $G4INSTALL/bin/geant4.sh"
    source $G4INSTALL/bin/geant4.sh
fi

echo "SIMRUN:: Run $DC_RUN Event $DC_EVENT Generator $CONFIG_RUN_TYPE Field $CONFIG_FIELD Energy $CONFIG_ENERGY Physicslist $CONFIG_PHYSICSLIST"

if [ -e "base_powheg.input" ]; then
  echo  "  >>> Run POWHEG <<<"
  chmod u+x runPowheg.sh
  nevts=1000
  if [ ! -z $DC_EVENT ]; then
    nevts=$DC_EVENT
  fi
  ./runPowheg.sh -v W $(expr $nevts \* 20) $CONFIG_SEED > powheg.log 2>&1
fi

runcommand "SIMULATION" "fastPOWHEG.C" sim.log 5
mv syswatch.log simwatch.log

exit 0

