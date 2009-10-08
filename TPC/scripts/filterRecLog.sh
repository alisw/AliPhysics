# marian.ivanov@cern.ch
# Filter error logs 
# Input   - errRec.logs  - text file with the list of error logs
# Output  - seg.log      - text file with the seg.faults runs
#         - segout.log   - text file with the output of the seg.fault runs
#                        - the err.log and out.log suppose to be in the same directory
# to get the list of logs:
# example:
# find /lustre/alice/alien/alice/data/2009/LHC09c_TPC/00008*/esdgsi3/*/err.log  > errRec.log

isOK=0
nonOK=0
rm seg.log
rm segout.log
for efile in `cat errRec.log`  ;do
 xxx=`cat $efile| grep segmentation`
 if [ -z "$xxx" ]
 then
  let isOK=isOK+1
  else
  let nonOK=nonOK+1
  echo nonOK=$nonOK
  echo "$efile" >>seg.log
  echo $efile
  ofile=`echo $efile| sed s_err_out_`
  cat $ofile >> segout.log
 fi     
done; 
#get the list
echo isOK=$isOK nonOK=$nonOK
