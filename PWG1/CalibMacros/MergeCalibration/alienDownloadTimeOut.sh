#
# sh script to download the data from alien to the local storage
#
# Arguments:
#   $1 fileList
#   $2 timeOut in seconds
#   $3 nprocess
#   $4 prefix
#
# Algorithm:
#   1. Check the input parameters
#   2. Group randomly the files into the group - generates nprocess download scripts
#   3. Start $nprocess to download the files - with forced timeout
# Example usage:
# 1. Create file list 
# alien_find /alice/data/2010/LHC10h/000138795/pass0_2 root_ar | grep root > file.list 
# 2. run script with timeout 60 and n processes 10 and the sarting diectory  -current directory
# $ALICE_ROOT/PWG1/CalibMacros/MergeCalibration/alienDownloadTimeOut.sh file.list 60 30 `pwd`
#

fileList=$1
timeOut=$2
nprocess=$3
prefix=$4
#
# 1. Check the input variables
#
if ! [ -e $fileList ]; then
  echo Input file does not exist
  exit;
fi; 

if [  ${#timeOut} -lt 1 ]; then
  echo Time out not specified
  exit;
fi; 

if [  ${#nprocess} -lt 1 ]; then
  echo N processes  not specified
  exit;
fi; 

#
#   2. Group randomly the files into the group - generates nprocess download scripts
#
nfiles=`wc $fileList | gawk '{print $1;}'`
echo NFiles"  "$nfiles"   "$fileList"  "TimeOut:$timeOut 
counterProcess=0;
mkdirhier tmpDownload
rm -rf tmpDownload/*.sh
#
# generate download scripts
#
for afile in `sort -r $fileList`; do
  dname=$prefix/`dirname $afile`
  fname=`basename $afile`
  if ! [ -e $dname/$fname ] ; then 
    mkdirhier $dname        
    echo "echo Date `date` >> job$counterProcess.log>> job$counterProcess.log  " >> tmpDownload/download$counterProcess.sh
    echo "echo alien_cp alien:////$afile  $dname 2>> job$counterProcess.log>> job$counterProcess.log  " >> tmpDownload/download$counterProcess.sh
    echo "timeout  $timeOut  alien_cp -t $timeOut -i alien.txt -m alien:////$afile  $dname 2>> job$counterProcess.log>> job$counterProcess.log  " >> tmpDownload/download$counterProcess.sh
    echo "cat alien.txt >> job$counterProcess.log  " >> tmpDownload/download$counterProcess.sh
    let counterProcess=$counterProcess+1
    if [ $counterProcess -gt $nprocess ]; then
      let counterProcess=0;
    fi; 
  fi;
done;

echo NFiles"  "$nfiles"   "$fileList"  "TimeOut:$timeOut 

#
#   3. Start $nprocess to download the files - with forced timeout
#
for job in `ls tmpDownload/download*.sh`; do 
  chmod u+x $job;
  echo $job
  command $job &
done; > download.log 
