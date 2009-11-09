#!/bin/sh
# 1 argument      - qname

# marian.ivanov@cern.ch
# filter corrupted files
# batch farm internaly used
#
# Assumption - the data organize in the Workspace - predefined directory structure
#            - each directory contain the list of files 
#            - raw.txt and esd.txt
# Output:
# esd.txt => esd.txt.Good and esd.txt.Bad 
# eaw.txt -> raw.txt.Good and raw.txt.Bad

qname=$1
mydir=`pwd`
for adir in `cat run.list`; do
    cd $adir
    up=`cat  raw.txt | grep -c .root`
    if [ $up -gt 0 ] ; then
	echo bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterRAW.C 
	bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterRAW.C
    fi;	
    cd $mydir
done     

mydir=`pwd`
for adir in `cat run.list`; do
    cd $adir
    up=`cat  esd.txt | grep -c .root`
    if [ $up -gt 1 ] ; then
	echo bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterESD.C
	bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterESD.C
    fi;
    cd $mydir
done
 

