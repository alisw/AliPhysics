#!/bin/bash

if [ -z "$1" ]
then
    echo "Verifies if the LDAP entry corresponds to the entry in the file."
    echo "Usage: verify.sh <3 DIGIT DET CODE>"
    echo "       verify.sh sys-[HLT|DCS|DAQ]"
    echo "       verify.sh global"
    echo "       verify.sh ALL"
    exit
fi

case "$1" in
    "global"  ) SEARCH="name=globalConfig"; FILE="Global.ldif";;
    "sys-HLT" ) SEARCH="system=HLT"; FILE="HLTsys.ldif";;
    "sys-DAQ" ) SEARCH="system=DAQ"; FILE="DAQsys.ldif";;
    "sys-DCS" ) SEARCH="system=DCS"; FILE="DCSsys.ldif";;
    "ALL"     ) LIST="EMC FMD HMP MCH MTR PHS PMD SPD SDD SSD TOF TPC TRD T00 V00 ZDC GRP ACO HLT global sys-HLT sys-DCS sys-DAQ"
                for I in $LIST
		do 
		  echo "Checking $I"
		  ./verify.sh $I
		done
		exit
		;;
    *         ) SEARCH="det=$1"; FILE="$1.ldif";;
esac

ldapsearch -H ldap://pcalishuttle02.cern.ch  -x -b "$SEARCH,o=shuttle_prod,dc=cern,dc=ch" -L -L -L > verify.out
cat verify.out | ./unfoldlines.pl > verify1.out

diff --ignore-space-change --ignore-blank-lines verify1.out $FILE | grep -v "> #" | grep -v '^0[a-z]'

rm verify.out verify1.out
