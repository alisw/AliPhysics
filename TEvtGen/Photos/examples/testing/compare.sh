#!/bin/bash   
# Written by N. Davidson (2008)
# This bash script demonstrates how to perform the "Analysis" 
# step from a directory outside MC-TESTER.
if test -z ${MCTESTERLOCATION}; then
	echo ""
	echo "ERROR: \$MCTESTERLOCATION not set. Analysis cannot be performed."
	echo "       Use: export MCTESTERLOCATION=<path_to_mc_tester>"
	echo ""
	exit -1
fi
#----------------------------------
# Change these variables
#----------------------------------
FILE1=$1 #first generation step file
FILE2=$2 #second generation step file
MCTESTER_DIR=${MCTESTERLOCATION} #location of MC-TESTER

#----------------------------------

MCTESTER_ANALYZE_DIR=${MCTESTER_DIR}/analyze
export MC_TESTER_LIBS_DIR=${MCTESTER_DIR}/lib

WORKING_DIR=`pwd`

#change to MCTester directory and run macros
cd $MCTESTER_ANALYZE_DIR 
root -b -q "ANALYZE.C(\"${WORKING_DIR}\",\"${WORKING_DIR}/${FILE1}\",\"${WORKING_DIR}/${FILE2}\")" 
root -b -q "BOOKLET.C(\"${WORKING_DIR}\")"
cd $WORKING_DIR

#copy base .tex file needed for booklet and create
cp ${MCTESTER_ANALYZE_DIR}/tester.tex ./ 
latex tester.tex
latex tester.tex 
latex tester.tex 
dvipdf tester

#do a bit of clean up
rm -rf tester.aux tester.log texput.log tester.toc 
rm -rf mc-results.aux booklet.aux 
