#!/bin/bash
export GCLIENT_SERVER_LIST="pcapiserv04.cern.ch:10000|pcapiserv05.cern.ch:10000|pcapiserv06.cern.ch:10000|pcapiserv07.cern.ch:10000"
echo ===========================
echo $PATH 
echo $ROOTSYS
echo $LD_LIBRARY_PATH
echo ==========================

root -b -x runBatch.C;

