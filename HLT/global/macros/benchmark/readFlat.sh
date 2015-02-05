#! /bin/bash
if [ $# -lt 1 ]
    then
			verbose=0
    else
				verbose=$1
    fi
if [ $# -lt 2 ]
    then
			filename="outFlatESD.dat"
    else
				filename=$2
    fi
aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/ReadFlatESD.C++'("'${filename}'",'${verbose}')' 2>&1| tee readFlat.out
