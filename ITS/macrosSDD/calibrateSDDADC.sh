#!/bin/sh
export period="LHC17k"
export filnam="CalibObjects"
export ocdbFile="Run0_999999999_v2_s0"
export dir="LHC17k"

for nrun in \
274690\

do

echo $nrun

root -l <<EOF
.x MakeSDDADCCalib.C+ ($nrun,"$dir","$filnam",2017,"$period",kTRUE);
.q
EOF

root -l <<EOF
.x MakeOCDBCorrectionFromOutputADCCalib.C+ ($nrun,"$dir","$period","$ocdbFile")
.q
EOF

done