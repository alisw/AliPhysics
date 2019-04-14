#!/bin/bash
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH

sys=$1
multi_start=$2
multi_end=$3
inputoptions=$4
optionnumber=$5


echo "========================================="
echo "############## PATH : ##############"
echo $PATH
echo "############## LD_LIBRARY_PATH : ##############"
echo $LD_LIBRARY_PATH
echo "############## ROOTSYS : ##############"
echo $ROOTSYS
echo "############## which root : ##############"
which root
echo "############## ALICE_ROOT : ##############"
echo $ALICE_ROOT
echo "############## which aliroot : ##############"
which aliroot
echo "############## system limits : ##############"
ulimit -a
echo "############## memory : ##############"
free 2> /dev/null || { [[ `uname` == Darwin ]] && top -l 1 -s 0 | head -8 | tail -3; }
echo "========================================="

root -b -q -x DrawXi1530.C+\($sys\,$multi_start\,$multi_end\,\"$inputoptions\",$optionnumber\) 2>&1
echo "############## memory after: ##############"
free -m


DATETIME=$(eval date)
TOKEN=582112631:AAEvRb36PE2lJoplyLoXK3NWAavezCX2f1Q
CHAT_ID=559725700
MESSAGE="$DATETIME Job Finished. Options: $1 [$2-$3] $4"
URL="https://api.telegram.org/bot$TOKEN/sendMessage"

#OUTPUT="$(ls AnalysisResults_Extracted*)"
#alien_cp -n $OUTPUT alien:///alice/cern.ch/user/b/blim/Xi1530Extractor/out/$OUTPUT

#curl -s -X POST $URL -d chat_id=$CHAT_ID -d text="$MESSAGE"