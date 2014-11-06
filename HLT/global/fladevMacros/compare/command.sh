#! /bin/bash
if [ $# -lt 1 ]
    then
	file1="flat.dat"
    else
	file1=$1
    fi
if [ $# -lt 2 ]
    then
	file2="normal.dat"
    else
	file2=$2
    fi
if [ $# -lt 3 ]
    then
	verbose="kFALSE"
    else
	verbose=$3
    fi



aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/CompareFlatESDs.C++'("'${file1}'","'${file2}'",'${verbose}')' | tee compare.out
