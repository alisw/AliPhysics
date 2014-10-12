#!/bin/bash

# set job and simulation variables as :
# ./Run.sh  --run <x> 

#
# Function to run a command and check return code 
function runcommand()
{
    echo -e "\n"
    echo -e "\n" >&2

    local type=$1 
    local scr=$2
    local log=$3 

    echo "* $type : $scr"
    echo "* $type : $scr" >&2
    
    time aliroot -b -q -x $scr 2> /dev/stdout | tee -a $log 
    local ext=$?
    local exp=${5-0}

    if test -f syswatch.log ; then 
	mv syswatch.log `basename $log .log`watch.log 
    fi 

    if [ "$ext" -ne "$exp" ]; then
        echo "*! $scr failed with exitcode $ext, expecting $exp"
        echo "*! $scr failed with exitcode $ext, expecting $exp" >&2
        echo "$scr failed with exitcode $ext, expecting $exp" \
	    > validation_error.message
        exit ${4-$ext}
    else
        echo "* $scr finished with the expected exit code ($exp), moving on"
        echo "* $scr finished with the expected exit code ($exp), moving on" >&2
    fi
}

#
# Function to clean up rec-points 
function cleanRecPoints()
{
    local alsoITS=${1-1}
    if test $alsoITS -gt 0 ; then 
	rm -f *.RecPoints.root
	return;
    fi
    for i in *.RecPoints.root ; do
	if test ! -f $i ; then continue ; fi 
	case $i in 
	    ITS*) ;; 
	    *) rm -f $i ;; 
	esac
    done
}
    
# Define the pt hard bin arrays
pthardbin_loweredges=( 0 5 11 21 36 57 84 117 152 191 234 )
pthardbin_higheredges=( 5 11 21 36 57 84 117 152 191 234 -1)

CONFIG_SEED=""
CONFIG_RUN_TYPE=""
CONFIG_BMIN=""
CONFIG_BMAX=""
CONFIG_QUENCHING=""
DC_RUN=""
DC_EVENT="1"
number=0
runAODTrain=0
runQATrain=0

while test "x$1" != "x"; do
    option="$1"
    shift
    
    case $option in 
        --run)		DC_RUN="$1";		shift ;;
        --event)	DC_EVENT="$1";		shift ;;
        --process)	CONFIG_RUN_TYPE="$1";	shift ;;
        --field)	: ; 			shift ;; # No-op
        --energy)	: ; 			shift ;; # No-op
        --physicslist)	: ;			shift ;; # No-op
        --bmin)		CONFIG_BMIN="$1";	shift ;;
        --bmax)		CONFIG_BMAX="$1";	shift ;;
        --pthardbin)	: ;			shift ;; # No-op
        --quench)	CONFIG_QUENCHING="$1";	shift ;;
	--sdd) 		: 			;; 	 # No-op
	--number)       number="$1";		shift ;;
	--qa)           runQATrain=1		      ;;
	--aod)          runAODTrain=1		      ;;
	*) echo "Unkown option: $option" >&2
    esac
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

mkdir -p input
test -f galice.root && mv galice.root ./input/galice.root
test -f Kinematics && mv Kinematics.root ./input/Kinematics.root
ls input

export CONFIG_SEED \
    CONFIG_RUN_TYPE \
    CONFIG_BMIN \
    CONFIG_BMAX \
    CONFIG_QUENCHING \
    DC_RUN \
    DC_EVENT

export ALIMDC_RAWDB1="./mdc1"
export ALIMDC_RAWDB2="./mdc2"
export ALIMDC_TAGDB="./mdc1/tag"
export ALIMDC_RUNDB="./mdc1/meta"

if [ -f "$G4INSTALL/bin/geant4.sh" ]; then
    echo "* Sourcing G4 environment from $G4INSTALL/bin/geant4.sh"
    source $G4INSTALL/bin/geant4.sh
fi

cat <<EOF

SIMRUN: Setup
  Seed:             $CONFIG_SEED
  Run:              $DC_RUN
  Events:           $DC_EVENT
  Event generator:  $CONFIG_RUN_TYPE
  Serial #:         $number
  Impact parameter: ${CONFIG_BMIN}-${CONFIG_BMAX}
  Quenching:        $CONFIG_QUENCHING
  Run QA:           $runQATrain
  Run AOD           $runAODTrain

Content of current directory:
EOF
ls -l 

echo "SIMRUN: Now read to process"

runcommand "SIMULATION"     "Simulate.C($DC_EVENT,$DC_RUN)" 	sim.log     5
runcommand "RECONSTRUCTION" "Reconstruct.C($DC_RUN)" 		rec.log    10
runcommand "TAG"            "Tag.C"         			tag.log    50
runcommand "CHECK" 	    "Check.C"    			check.log  60
if test $runQATrain -gt 0 ; then 
    runcommand "QA"         "QA.C($DC_RUN)"  		        qa.log    100
fi
if test $runAODTrain -gt 0 ; then 
    runcommand "AOD"        "AOD.C($DC_RUN)"  		        aod.log   100
fi

rm -f *.Hits.root *.Digits.root *.SDigits.root
cleanRecPoints 0

exit 0
#
# EOF
#

