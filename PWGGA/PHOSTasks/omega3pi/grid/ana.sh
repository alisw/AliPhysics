#!/bin/bash
export GCLIENT_SERVER_LIST="pcapiserv04.cern.ch:10000|pcapiserv05.cern.ch:10000|pcapiserv06.cern.ch:10000|pcapiserv07.cern.ch:10000"
export GRID_TOKEN=OK
echo ===========================
echo $PATH 
echo $ROOTSYS
echo $LD_LIBRARY_PATH
echo $GRID_TOKEN
echo dataset = $1.xml
echo ==========================

aliroot -b -q AnaTaskOmega3pi.C\(\"$1.xml\"\)
echo "########################    Train finished    ###########################"

