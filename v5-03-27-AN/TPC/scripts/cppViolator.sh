#!/bin/bash
#
# simple script to filter  warnings in TPC code
# root warning are filtered out
#
# Usage: 
#    cd $ALICE_ROOT
#    $ALICE_ROOT/TPC/scripts/cppViolator.sh
# warnings.txt file will be created 
# Please attach the the warnings.txt to the request for commit
#


make clean-TPC
make -i -j 5 &> compile.txt

cat compile.txt | grep -v $ROOTSYS | grep -v G__ | grep -v TMatrix | grep -v TMap|  grep warning > warnings.txt

cat warnings.txt





