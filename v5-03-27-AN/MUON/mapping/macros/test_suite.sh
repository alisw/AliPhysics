#!/bin/sh
# $Id$

# Test suite for running St12 mapping test macros in a batch mode.
# The macros outputs are saved in test_out directory.
#
# by I. Hrivnacova, IPN Orsay

OUTDIR="test_out"
CURDIR=`pwd`

cd $ALICE_ROOT/MUON/mapping/macros
if [ ! -d $OUTDIR ] 
then
  mkdir -p $OUTDIR 
fi

for TEST in `ls testSt12*C` `ls testSt345*C` testDE.C timeMapping.C
do
  TESTNAME=`echo $TEST | sed 's/\.C//g'`
  echo "Running $TESTNAME ..."
  #aliroot -b >& testSt12AllIndices.out << EOF
  aliroot -b >& $OUTDIR/$TESTNAME".out" << EOF
  .L $TEST+
  $TESTNAME();
  .q
EOF

done

mv testExistingPads.*.out $OUTDIR

cd $CURDIR
