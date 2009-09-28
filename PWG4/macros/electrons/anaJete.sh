#!/bin/bash

#
# Author: K. Read
#

 export anaInputData=ESD
 export INDIR=/work2/data/test
 export PATTERN=
 export NFILES=1

#export anaInputData=AOD
#export INDIR=/work2/data/LHC08d6/AOD/000
#export PATTERN=01
#export NFILES=1

#export anaInputData=AOD
#export INDIR=/work2/data/bjetfilter/AOD/117005
#export PATTERN=00
#export NFILES=5

 export MODE=0

 export CONFIG2=ConfigJetAnalysisFastJet.C
 export CONFIG3=ConfigAnalysisElectron

 export SIMPATH=/work2/data/bjetfilter/AOD/117005/004
#export SIMPATH=/work2/data/LHC09a5/000042
#export SIMPATH=/work2/data/LHC08d6/AOD/000/010

 export SEVENT=0

root -q -b -l anaJete.C
