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
        ocdb="local://$OCDB10"
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
  rm galice.root
  mkdir ${iFile}
  
  if [ $config = "flat" ]
	then
	  echo "using flat ESD converter"
	  aliroot -q -b ../config_Flat.C'("'${iFile}'","'${ouputfileName}'")' $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=RootWriter  ignore-hltout")' 2>&1|tee recraw-local.log -a
	  cat ${iFile}/* >> outFlatESD.dat
	else
	  echo "using normal ESD converter"
	  aliroot -b -q -l $ALICE_ROOT/HLT/exa/recraw-local.C'("'${file}'","'${ocdb}'", '${start}', '${end}', "HLT","chains=GLOBAL-esd-converter ignore-hltout")' 2>&1|tee recraw-local.log -a
	  mv AliESD* ${iFile}/
	fi
  
  iFile=$((iFile+1))
done
