#!/bin/bash
set -e

if [[ $# -lt 2 ]]; then
    echo "You must specify at least two parameters where to write output."
    exit -1
fi

EPOSFILE1=${1}
EPOSFILE2=${2}

echo "running in `pwd`, writing EPOSroot to ${EPOSFILE1} and ${EPOSFILE2}"

# prepare environment
source /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv EPOS::v3.111-8)

# prepare EPOS configuration
OptnsEpos3="pp"
ENERGY="13000"
MAXEPOSEVENTS=100
EPOSNFREEZE=10

cp ${ALICE_PHYSICS}/PWG/MCLEGO/EPOS/${OptnsEpos3}.optns .
sed -i "s/\(istmax\)\(\s\+[0-9]\+\)/\1 10/g" ${OptnsEpos3}.optns
sed -i "s/\(ecms\)\(\s\+[0-9]\+\)/\1 $ENERGY/g" ${OptnsEpos3}.optns
sed -i "s/\(nfull\)\(\s\+[0-9]\+\)/\1 $(( $MAXEPOSEVENTS / $EPOSNFREEZE ))/g" ${OptnsEpos3}.optns
sed -i "s/\(nfreeze\)\(\s\+[0-9]\+\)/\1 $EPOSNFREEZE/g" ${OptnsEpos3}.optns

# generate events and alternate between two files
until true; do
    epos $OptnsEpos3 > epos.out 2> epos.err;
    if [ ! -f $EPOSFILE1 ]; then
	mv z-${OptnsEpos3}.root $EPOSFILE1;
    elif [ ! -f $EPOSFILE2 ]; then
	mv z-${OptnsEpos3}.root $EPOSFILE2;
    else
	sleep 1s
    fi
done
