#! /bin/bash
if [ $# -lt 1 ]
    then
			filename="outFlatESD.dat"
    else
				filename=$1
    fi
aliroot -b -l -q $ALICE_ROOT/HLT/global/LoadLibs.C $ALICE_ROOT/HLT/global/ReadFlatESD.C++'("'${filename}'")' | tee readFlat.out
