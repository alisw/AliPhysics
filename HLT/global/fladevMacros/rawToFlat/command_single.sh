#! /bin/bash
if [ $# -lt 1 ]
    then
			inputFileName="raw.root"
    else
				inputFileName=$1
    fi
if [ $# -lt 2 ]
    then
			outputDirectory="outFlat"
    else
				outputDirectory=$2
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

aliroot -q -b ../config_Flat.C'("'${outputDirectory}'", "'${outputFileName}'" )' $ALICE_ROOT/HLT/exa/recraw-local.C'("'${inputFileName}'","'${ocdb}'", '${start}', '${end}', "HLT", "chains=RootWriter  ignore-hltout")' |tee output.out -a

printf "now merging files"
cat ${outputDirectory}/* > outFlatESD.dat
