#! /bin/bash
if [ $# -lt 1 ]
    then
			echo "please specify directory with input raw data"
			exit
    else
				fileList=$1
    fi
if [ $# -lt 2 ]
    then
			config="normal"
    else
				config=$2
    fi
if [ $# -lt 3 ]
    then
			ouputfileName="outFlatHLT.dat"
    else
				ouputfileName=$3
    fi
if [ $# -lt 4 ]
    then
        ocdb="local://$OCDB10local"
    else
				ocdb=$4
    fi
if [ $# -lt 5 ]
    then
        start=-1
    else
				start=$5
    fi
if [ $# -lt 6 ]
    then
        end=-1
    else
				end=$6
    fi
        
aliroot ../createHistos.C -q
iFile=1

for file in ${fileList}*/raw.root
do
  #  dir=${dir%*/}
  #  echo ${dir##*/}
  echo "Now processing ${file}"
  mkdir ${iFile}
  
  if [ $config = "flat" ]
	then
	  echo "using flat ESD converter"
	  aliroot -q -b ../config_Flat.C'("'${iFile}'","'${ouputfileName}'")' $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=RootWriter  ignore-hltout")' 2>&1|tee recraw-local.log -a
	  cat ${iFile}/*.dat >> outFlatESD.dat
	else
	  echo "using normal ESD converter"
	  aliroot -b -q -l $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT","chains=GLOBAL-esd-converter ignore-hltout")' 2>&1|tee recraw-local.log -a
	  mv AliESD* ${iFile}/
	fi
	mv syswatch.log ${iFile}/syswatch.log 
	if [ $iFile -eq 1 ]
		then
			cp 1/syswatch.log syswatch_merged.log
		else
			sed 1d $iFile/syswatch.log >> syswatch_merged.log
		fi
	mv *.root ${iFile}/
  iFile=$((iFile+1))
done
