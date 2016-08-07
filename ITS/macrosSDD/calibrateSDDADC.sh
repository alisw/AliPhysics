export nrun=246424
export period="15o"
export dir="LHC15o"
export filnam="QAresults_barrel_246424.root"
export ocdbFile="Run246393_999999999_v2_s0"

root -l <<EOF
.x MakeSDDADCCalib.C+ ($nrun,"$dir","$filnam");
.q
EOF

root -l <<EOF
.x MakeOCDBCorrectionFromOutputADCCalib.C+ ($nrun,"$dir","$period","$ocdbFile")
.q
EOF
