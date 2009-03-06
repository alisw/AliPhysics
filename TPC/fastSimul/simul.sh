#!/bin/sh

# 1 argument      - the path to the environment setup
# 2 argument      - the job ID
# 3 argument      - number of events in the file
# 4 argument      - output path

# Example
# myvar=0
# $ALICE_ROOT/TPC/fastSimul/simul.sh  /u/miranov/.balice64HEAD0108 $myvar 1000 `pwd`

#
# 1 SETUP given ROOT and ALIROOT
#
echo   $1
source $1
echo  $ROOTSYS
which root.exe
which aliroot
#
#  make directory
#

cd $4
mkdir $2
cd $2
cp ~/rootlogon.C .
echo Job ID  $2
echo
echo PWD `pwd`

command aliroot  -q -b  "$ALICE_ROOT/TPC/fastSimul/simul.C($3)"