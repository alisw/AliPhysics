#! /bin/bash
if [ $# -lt 1 ]
    then
	fileList=$PWD
    else
	fileList=$1
    fi
if [ $# -lt 2 ]
    then
	useFriends=1
    else
	useFriends=$2
    fi
if [ $# -lt 3 ]
    then
	useHLTtree=0
    else
	useHLTtree=$3
    fi
if [ $# -lt 4 ]
    then
	verbose=0
    else
	verbose=$4
    fi
    
    iFile=1
for file in ${fileList}*/AliESDs.root
do
  mkdir ${iFile}
  replace="AliESDfriends"
  fileFriends="${file/AliESDs/$replace}"
  aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/FlatESDConverter.C++'("'${file}'", "'${fileFriends}'", "'${iFile}'/out.dat",'${useFriends}', '${useHLTtree}','${verbose}')' 2>&1| tee convert.log -a
  cat ${iFile}/out.dat >> outFlatESD.dat
  mv syswatch.log ${iFile}/syswatch.log
	if [ $iFile -eq 1 ]
		then
			cp 1/syswatch.log syswatch_merged.log
		else
			sed 1d $iFile/syswatch.log >> syswatch_merged.log
	fi
  iFile=$((iFile+1))
done
