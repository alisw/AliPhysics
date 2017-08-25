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

export ARG="PSJetTask_merge.C(\"$1\",$2)"
root -b -q -x $ARG
RET=$?
if [ "$RET" != "0" ];then
  echo "======== ERROR : PSJetTask.C finished with NON zero code: $RET ========"
  if [ "$RET" -gt 128 ] && [ "$RET" -lt 160 ]; then
    let sig="$RET - 128"
    sigs='HUP INT QUIT ILL TRAP ABRT BUS FPE
    KILL USR1 SEGV USR2 PIPE ALRM TERM STKFLT
    CHLD CONT STOP TSTP TTIN TTOU URG XCPU
    XFSZ VTALRM PROF WINCH IO PWR SYS'
    sig=SIG`echo $sigs | awk '{ print $'"$sig"' }'`
    echo "======== it appears to have been killed with signal: $sig ========"
  fi
  exit $RET
fi

echo "======== PSJetTask_merge.C finished with exit code: $? ========"
echo "############## memory after: ##############"
free -m
