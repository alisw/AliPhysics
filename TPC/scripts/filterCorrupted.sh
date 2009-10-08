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
    bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterRAW.C
    cd $mydir
    done;     
done;
#
mydir=`pwd`
for adir in `cat run.list`; do
cd $adir
echo bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterESD.C
bsub -q $qname  aliroot -b -q  $ALICE_ROOT/TPC/macros/filterESD.C
cd $mydir
done;


