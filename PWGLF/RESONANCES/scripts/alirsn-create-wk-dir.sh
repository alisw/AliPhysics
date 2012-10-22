#!/bin/bash

MY_ALIRSN_DIR="$0"
MY_ALIRSN_DIR=`dirname $MY_ALIRSN_DIR`
MY_ALIRSN_DIR=`dirname $MY_ALIRSN_DIR`

FILES="AddAMRsn.C RunALICE.C SetupAnalysisPlugin.C"

for FILE in $FILES;do
  cp $MY_ALIRSN_DIR/PWGLF/RESONANCES/macros/lego_train/$FILE .
done