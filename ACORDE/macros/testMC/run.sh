#!/bin/sh
echo '***** Running simulation and reconstruction Test for ACORDE Detector *****'
echo '***** Comments to: Mario Rodriguez Cahuantzi *****'
echo '***** <mrodrigu@mail.cern.ch>, <mario.rocah@gmail.com> *****'
echo '***** Comments to: Pedro Podesta <podesta@cern.ch> *****'
echo '***** FCFM - BUAP, Puebla, Pue. Mexico 2010 *****'
echo '***** Current Date *****'
date
echo '************************'
echo '***** Cleaning *****'
rm -rf *.ps *.root *.dat *.log fort* hit hough raw* recraw/*.root recraw/*.log recraw/*.eps recraw/*.ps rm AliHL* raw* HLT GRP *.eps
echo '***** Running cosmic simulation *****'
aliroot -b -q  cosmicSim.C      2>&1 | tee simCosmic.log
echo '***** Running cosmic reconstruction *****'
aliroot -b -q  rec.C      2>&1 | tee recCosmic.log
echo '***** Running cosmic reconstruction with RAW Input *****'
cd recraw
ln -s ../raw.root
aliroot -b -q rec.C	  2>&1 | tee recCosmic.log
echo '***** Current Date *****'
date
echo '************************'
echo '***** Test for ACORDE Concluded *****'
echo '***** FCFM-BUAP, Puebla, Pue. Mexico *****'
echo '***** MRC: <mrodrigu@mail.cern.ch>, <mario.rocah@gmail.com> *****'
