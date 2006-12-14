#!/bin/bash
export GCLIENT_SERVER_LIST="pcapiserv01.cern.ch:10000|pcapiserv02.cern.ch:10000"
echo ===========================
echo $PATH 
echo $ROOTSYS
echo $LD_LIBRARY_PATH
echo ==========================

root -b -x runAnalysis.C;

