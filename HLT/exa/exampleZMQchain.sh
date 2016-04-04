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

histograms=( "hXTRKVtx" "hYTRKVtx" "hZTRKVtx" "hXTRKDefMult" "hYTRKDefMult"  "hZTRKDefMult" )
sourcePID=()
for hist in ${histograms[@]}; do
  echo "starting histogram source ${hist}"
  ZMQhistSource name=${hist} entries=10 sleep=0.01 distribution="exp(-0.5*((x-0.)/0.1)**2)" range="'-0.5,0.5'" out="PUSH>ipc:///tmp/example_proxyin" &>/dev/null &
  sourcePID+=(${!})
done

ZMQproxy  in="PULL@ipc:///tmp/example_proxyin" out="PUSH@ipc:///tmp/example_proxyout"&
sourcePID+=(${!})
ZMQROOTmerger in="PULL>ipc:///tmp/example_proxyout" out="PUB@ipc:///tmp/example_mergerout" MaxObjects=10 pushback-period=1&
sourcePID+=(${!})
ZMQhistViewer in="SUB>ipc:///tmp/example_mergerout" &
sourcePID+=(${!})

sleep 1

while true; do sleep 1; done
cleanupAndExit
