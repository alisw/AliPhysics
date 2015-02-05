#! /bin/bash
if [ $# -lt 1 ]
    then
			echo "please specify directory with input raw data"
			exit
    else
				fileList=$1
    fi
config=${2:-"c"}
ocdb=${3:-"local://$OCDB10"}
start=${4:--1}    
end=${5:--1}
outpufileName=${6:-"outFlatHLT.dat"}        
iFile=1

for file in ${fileList}*/raw.root
do
  #  dir=${dir%*/}
  #  echo ${dir##*/}
  echo "Now processing ${file}"
  mkdir ${iFile}
  
  if [ $config = "rawtoflat" -o $config = "f" ]; then
	  echo "running raw->flatESD conversion"
	  aliroot -q -b ../config.C'("GLOBAL-flat-esd-converter", "'${iFile}'","'${ouputfileName}'")' $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=RootWriter  ignore-hltout")' 2>&1|tee recraw-local.log -a
	  cat ${iFile}/*.dat >> outFlatESD.dat
  elif [ $config = "rawtoesdtoflat" -o $config = "ef" ]; then
	  echo "running raw->Esd->flatESD conversion"
	  aliroot -q -b ../config.C'("esd-to-flat-conversion", "'${iFile}'","'${ouputfileName}'")' $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=RootWriter  ignore-hltout")' 2>&1|tee recraw-local.log -a
	  cat ${iFile}/*.dat >> outFlatESD.dat
	elif [ $config = "compare" -o $config = "c" ]; then
	  echo "running raw->flatESD AND raw->Esd->flatESD conversion, compare outputs"
	  aliroot -q -b $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=compare-flat  ignore-hltout")' 2>&1|tee recraw-local.log -a
	elif [ $config = "rawtoesd" -o $config = "e" ]; then
	  echo "running raw->ESD conversion"
	  aliroot -b -q -l $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT","chains=GLOBAL-esd-converter ignore-hltout")' 2>&1|tee recraw-local.log -a
	  mv AliESD* ${iFile}/
	else
			echo "invalid option $(config)"  
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



