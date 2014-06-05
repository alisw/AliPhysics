#! /bin/bash
if [ $# -lt 1 ]
    then
			echo "please specify directory with input ESD data"
			exit
    else
				fileList=$1
    fi
    
    iFile=1
for file in ${fileList}*/AliESDs.root
do
  replace="AliESDfriends"
  fileFriends="${file/AliESDs/$replace}"
  aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/FlatESDConverter.C++'("'${file}'", "'${fileFriends}'", "'${iFile}'.dat",kFALSE, kTRUE)' 2>&1| tee convert.log
  cat ${iFile}.dat >> outFlatESD.dat
  
  iFile=$((iFile+1))
done
