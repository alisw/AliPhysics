#/bin/sh
#set -x

# Script for processing DETdigits.C (or digitITSDET.C) macros 
# and generating files with histograms.
# Generated histograms can be then plot with plotDigits2.C.
# Currently the macro has to be run in the directory
# with the aliroot output, and one has to specify the 
# number of events and the number of files in the chain
# to be prcessed.
#
# By E.Sicking, CERN; I. Hrivnacova, IPN Orsay


scriptdir=$ALICE_ROOT/test/vmctest/scripts

export NEVENTS=10
export NFILES=4

for DET in SDD SSD SPD TPC TRD TOF EMCAL HMPID PHOS; do
  echo Processing $DET digits
  aliroot -b -q $scriptdir/digits${DET}.C\($NEVENTS,$NFILES\)
done

