#!/bin/bash
export GCLIENT_SERVER_LIST="pcapiserv04.cern.ch:10000|pcapiserv05.cern.ch:10000|pcapiserv06.cern.ch:10000|pcapiserv07.cern.ch:10000"
export GRID_TOKEN=OK
echo ===========================
echo $PATH 
echo $ROOTSYS
echo $LD_LIBRARY_PATH
echo $GRID_TOKEN
echo ==========================

aliroot -b -q AnaPi0Select.C
echo "########################    Train finished    ###########################"

