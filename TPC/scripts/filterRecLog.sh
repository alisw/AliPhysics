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
rm abort.log
rm abortout.log
rm syswatchAbort.log
rm syswatchAll.log
echo hname/C:sname/C:id0/I:id1/I:id2/I:first/I:stampSec/I:mi.fMemUsed/F:mi.fSwapUsed/F:cI.fUser/F:cI.fSys/F:pI.fMemResident/F:pI.fMemVirtual/F:pI.fCpuUser/F:pI.fCpuSys/F:stampOldSec/I:miOld.fMemUsed/F:miOld.fSwapUsed/F:cIOld.fUser/F:cIOld.fSys/F:pIOld.fMemResident/F:pIOld.fMemVirtual/F:pIOld.fCpuUser/F:pIOld.fCpuSys/F > syswatchAbort.log

echo hname/C:sname/C:id0/I:id1/I:id2/I:first/I:stampSec/I:mi.fMemUsed/F:mi.fSwapUsed/F:cI.fUser/F:cI.fSys/F:pI.fMemResident/F:pI.fMemVirtual/F:pI.fCpuUser/F:pI.fCpuSys/F:stampOldSec/I:miOld.fMemUsed/F:miOld.fSwapUsed/F:cIOld.fUserF:cIOld.fSys/F:pIOld.fMemResident/F:pIOld.fMemVirtual/F:pIOld.fCpuUser/F:pIOld.fCpuSys/F > syswatchAll.log

#
for efile in `cat errRec.log`  ;do
 xxx=`cat $efile| grep segmentation`
 xxx=$xxx`cat $efile| grep Aborted`
 sysfile=`echo $efile| sed s_err.log_syswatch.log_`
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
  cat $ofile   >> abortout.log
  cat $sysfile | grep -v hname\/C:sname\/C: >> syswatchAbort.log
 fi
  cat $sysfile | grep -v hname\/C:sname\/ >> syswatchAll.log
done; 

#get the list
echo isOK=$isOK nonOK=$nonOK
