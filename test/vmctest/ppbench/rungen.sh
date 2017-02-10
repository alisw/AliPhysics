#!/bin/sh
# $Id$

NEVENTS=100
OUTDIR=gen

rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
aliroot -b -q rungen.C\($NEVENTS\)      2>&1 | tee gen.log
rm -fr $OUTDIR
mkdir $OUTDIR
mv galice.root Kinematics*.root gen.log $OUTDIR

rm -rf *.root *.dat *.log fort* hlt hough raw* recraw/*.root recraw/*.log
