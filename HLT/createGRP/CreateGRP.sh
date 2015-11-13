#!/bin/bash
export DIM_DNS_HOST=aldaqecs.cern.ch
export DIM_DNS_NODE=alidcsdimdns.cern.ch

export ALIHLT_T_HCDBDIR=./HCDB

RUN_NUMBER=169838
DETECTOR_LIST=TPC,SPD,SSD
BEAM_TYPE=pp
RUN_TYPE=PHYSICS
#These values for START_SHIFT and END_SHIFT are copied from PendolinoDriver.sh
let START_SHIFT=`date +%s`-300
END_SHIFT=`date +%s`

USE_DEFAULT_VALUES=kTRUE

echo Starting ALICE GRP Creation with Run Number $RUN_NUMBER, Detector List $DETECTOR_LIST, Beam Type $BEAM_TYPE, Run Type $RUN_TYPE

aliroot -l -q -b "AliHLTCreateGRP.C($RUN_NUMBER, \"$DETECTOR_LIST\", \"$BEAM_TYPE\", \"$RUN_TYPE\", $START_SHIFT, $END_SHIFT, $USE_DEFAULT_VALUES)"

retVal=$?

if [ $retVal -eq 0 ]; then
    echo "GRP Creation finished successfully"
else
    echo "GRP Creation failed"
fi
