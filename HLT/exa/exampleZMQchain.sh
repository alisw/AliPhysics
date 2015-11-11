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

sourceout=${1-"PUSH>ipc:///tmp/example_mergerin"}

histograms=( "hXTRKVtx" "hYTRKVtx" "hZTRKVtx" "hXTRKDefMult" "hYTRKDefMult"  "hZTRKDefMult" )
sourcePID=()
for hist in ${histograms[@]}; do
  echo "starting histogram source ${hist}"
  ZMQhistSource name=${hist} entries=10 sleep=0.01 distribution="exp(-0.5*((x-0.)/0.1)**2)" range="'-0.5,0.5'" out=${sourceout} &>/dev/null &
  sourcePID+=(${!})
done

ZMQROOTmerger in="PULL@ipc:///tmp/example_mergerin" out="PUB@ipc:///tmp/example_mergerout" MaxObjects=10 pushback-period=100&
sourcePID+=(${!})
ZMQhistViewer in="SUB>ipc:///tmp/example_mergerout" &
sourcePID+=(${!})

sleep 1

while true; do sleep 1; done
cleanupAndExit
