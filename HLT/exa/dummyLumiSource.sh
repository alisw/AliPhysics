#!/usr/bin/env bash
cleanupAndExit()
{
  for pid in ${sourcePID[@]}; do
    echo "killing $pid"
    kill ${pid}
  done
  exit
}

trap cleanupAndExit SIGINT
trap cleanupAndExit SIGTERM
trap cleanupAndExit SIGQUIT
trap cleanupAndExit SIGHUP

[[ -z $ALIHLT_CORE_DIR ]] && echo "ALIHLT_CORE_DIR not defined!" && exit
sourcePath=$ALIHLT_CORE_DIR/monitoring/ZMQproxy/bin
out=${1-"PUSH>ipc:///tmp/lumimergerin"}
monitor=${2-"SUB>ipc:///tmp/lumimergermonitor"}
sourcePID=()

ZMQhistSource name="hXTRKVtx" entries=10 sleep=0.001     distribution="exp(-0.5*((x+0.3)/0.02)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})
ZMQhistSource name="hYTRKVtx" entries=10 sleep=0.001     distribution="exp(-0.5*((x+0.2)/0.02)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})
ZMQhistSource name="hZTRKVtx" entries=10 sleep=0.001     distribution="exp(-0.5*((x+0.1)/0.02)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})
ZMQhistSource name="hXTRKDefMult" entries=10 sleep=0.001 distribution="exp(-0.5*((x-0.1)/0.06)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})
ZMQhistSource name="hYTRKDefMult" entries=10 sleep=0.001 distribution="exp(-0.5*((x-0.2)/0.06)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})
ZMQhistSource name="hZTRKDefMult" entries=10 sleep=0.001 distribution="exp(-0.5*((x-0.3)/0.06)**2)" range="'-0.7,0.7'" out=${out} &>/dev/null &
sourcePID+=(${!})

ZMQROOTmerger in="PULL@ipc:///tmp/lumimergerin" out="PUB@ipc:///tmp/lumimergerout" mon="REP@ipc:///tmp/lumimergermon" -pushback-period=5000 &
sourcePID+=(${!})

ZMQDIMlumiregServer -s "REQ>ipc:///tmp/lumimergerout" -publishperiod 5000 -dimdns aldaqecs.cern.ch -vertexservicename VERTEX -servername ALICE_HLT &
sourcePID+=(${!})

while true; do sleep 1; done
cleanupAndExit
