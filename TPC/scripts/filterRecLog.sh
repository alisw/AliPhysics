# marian.ivanov@cern.ch
# Filter error logs 
# Input   - errRec.logs  - text file with the list of error logs
# Output  - abort.log      - text file with the seg.faults runs
#         - abortout.log   - text file with the output of the seg.fault runs
#                          - the err.log and out.log suppose to be in the same directory
#         -syswatchAbort.log
# 
isOK=0 
nonOK=0
#
rm -f abort.log
rm -f abortout.log
rm -f syswatchAbort.log
rm -f syswatchAll.log
rm -f seg0.out

echo hname/C:sname/C:id0/I:id1/I:id2/I:first/D:stampSec/D:mi.fMemUsed/D:mi.fSwapUsed/D:cI.fUser/D:cI.fSys/D:pI.fMemResident/D:pI.fMemVirtual/D:pI.fCpuUser/D:pI.fCpuSys/D:stampOldSec/D:miOld.fMemUsed/D:miOld.fSwapUsed/D:cIOld.fUser/D:cIOld.fSys/D:pIOld.fMemResident/D:pIOld.fMemVirtual/D:pIOld.fCpuUser/D:pIOld.fCpuSys/D > syswatchAbort.log
echo hname/C:sname/C:id0/I:id1/I:id2/I:first/D:stampSec/D:mi.fMemUsed/D:mi.fSwapUsed/D:cI.fUser/D:cI.fSys/D:pI.fMemResident/D:pI.fMemVirtual/D:pI.fCpuUser/D:pI.fCpuSys/D:stampOldSec/D:miOld.fMemUsed/D:miOld.fSwapUsed/D:cIOld.fUser/D:cIOld.fSys/D:pIOld.fMemResident/D:pIOld.fMemVirtual/D:pIOld.fCpuUser/D:pIOld.fCpuSys/D > syswatchAll.log


#
for efile in `cat errRec.log`  ;do
 xxx=`cat $efile| grep segmentation`
 xxx=$xxx`cat $efile| grep Aborted`
 xxx=$xxx`cat $efile| grep floating`
 sysfile=`echo $efile| sed s_err_syswatch_`
 # 
 if [ -z "$xxx" ]
 then
  let isOK=isOK+1
  else
  let nonOK=nonOK+1
  echo nonOK=$nonOK
  echo "$efile" >>abort.log
  echo $efile
  ofile=`echo $efile| sed s_err_out_`
  cat $efile   >> abortout.log
  cat $ofile   >> abortout.log
  #cat $sysfile | grep -v hname\/C:sname\/C: >> syswatchAbort.log
 fi
  #cat $sysfile | grep -v hname\/C:sname\/ >> syswatchAll.log 
done; 

#
# connection problems
#
rm networkProblem.log
touch networkProblem.log
netOK=0;
netNonOK=0;
for efile in `cat errRec.log`  ;do
    xxx=`cat $efile| grep tcp_connect`
    if [ -z "$xxx" ]
    then
      let netOK=netOK+1
    else
      let netNonOK=netNonOK+1
      echo $efile >> networkProblem.log 
    fi;
done;

rm nfsProblem.log
touch nfsProblem.log
nfsOK=0;
nfsNonOK=0;
for efile in `cat errRec.log`  ;do
    xxx=`cat $efile| grep tcp_connect`
    xxx=$xxx`cat $efile| grep  Stale\ NFS\ file\ handle`
    if [ -z "$xxx" ]
    then
      let nfsOK=nfsOK+1
    else
      let nfsNonOK=nfsNonOK+1
      echo $efile >> nfsworkProblem.log 
    fi;
done;

#
# Print stat
#
echo isOK=$isOK nonOK=$nonOK
echo netOK=$netOK netNonOK=$netNonOK
echo nfsOK=$nfsOK netNonOK=$nfsNonOK
#
# filter segmentation fault
#
rm seg0.out
for a in `cat  abort.log |sed s_err_out_ ` ;do
    cat $a | grep 0x | grep \# >> seg0.out
done;

for a in `cat  abort.log ` ;do
    cat $a | grep 0x | grep \# >> seg0.out
done;
