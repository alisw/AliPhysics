#!/bin/bash

#
# Author: K. Read
#

export INDIR=/work2/data/test
export PATTERN=
#export INDIR=/work2/data/LHC08d6/AOD/000
#export PATTERN=01
export NFILES=1
export MODE=0
export CONFIG2=ConfigJetAnalysisFastJet.C
export CONFIG3=ConfigAnalysisElectron
export SIMPATH=/work2/data/bjetfilter/AOD/117005/004
export SEVENT=0

root -q -b -l anaJete.C
