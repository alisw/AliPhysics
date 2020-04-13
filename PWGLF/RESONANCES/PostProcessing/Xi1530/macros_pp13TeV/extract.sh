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