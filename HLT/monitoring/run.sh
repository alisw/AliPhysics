#!/bin/bash

#tpc a00
POLLNODE=cn020.internal
POLLPORT=49153
SERVERPORT=49001
CDBPATH="local://$ALIHLT_HCDBDIR"
POLLPERIOD=0.1
#compression monitoring
CONFIG=2

RUN_DIR=/tmp/mon_chain
mkdir -p $RUN_DIR
pushd $RUN_DIR

MACRO_DIR=$(dirname $(readlink -f $0))
ln -sf ${MACRO_DIR}/rootlogon.C .
ln -sf ${MACRO_DIR}/MonitorSandbox.C .
run_compression_monitor(){
    aliroot -b -q -l rootlogon.C MonitorSandbox.C+"(\"${POLLNODE}\",${POLLPORT},${SERVERPORT},${POLLPERIOD},\"${CDBPATH}\",${CONFIG})"
}

while true; do
    run_compression_monitor
    sleep 10
done | tee mon_chain.log

popd
