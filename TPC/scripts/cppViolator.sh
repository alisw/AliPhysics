#!/bin/bash
#
# simple script to filter  warnings in TPC code
# root warning are filtered out
# 


make clean-TPC
make -i -j 5 &> compile.txt

cat compile.txt | grep -v $ROOTSYS | grep warning




