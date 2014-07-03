#! /bin/bash
if [ $# -lt 1 ]
    then
			echo "please specify directory with input raw data"
			exit
    else
				fileList=$1
    fi
iFile=1

for file in ${fileList}*/raw.root
do
  #  dir=${dir%*/}
  #  echo ${dir##*/}
  echo "Now processing ${file}"
  mkdir ${iFile}

	 echo "using normal ESD converter"
	 aliroot -b -q ../rec.C'("'${file}'")' 2>&1|tee rec.log -a
	 mv *.root ${iFile}/


	mv syswatch.log ${iFile}/syswatch.log 
	if [ $iFile -eq 1 ]
		then
			cp 1/syswatch.log syswatch_merged.log
		else
			sed 1d $iFile/syswatch.log >> syswatch_merged.log
		fi
  iFile=$((iFile+1))
done
