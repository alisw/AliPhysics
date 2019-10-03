#!/bin/bash
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
echo "========================================="
echo "############## PATH : ##############"
echo $PATH
echo "############## LD_LIBRARY_PATH : ##############"
echo $LD_LIBRARY_PATH
echo "############## ROOTSYS : ##############"
echo $ROOTSYS
echo "############## which root : ##############"
which root
echo "############## ALICE_ROOT : ##############"
echo $ALICE_ROOT
echo "############## which aliroot : ##############"
which aliroot
echo "############## system limits : ##############"
ulimit -a
echo "############## memory : ##############"
free -m
echo "========================================="

aliroot -b << EOF
gSystem->SetIncludePath("-I$ALICE_ROOT/MUON");
.L QAMerge.C+
QAMerge("$1","VAR_MERGED_OUTPUT_NAME");
EOF

echo "======== QAmerge.C finished with exit code: $? ========"
echo "############## memory after: ##############"
free -m
