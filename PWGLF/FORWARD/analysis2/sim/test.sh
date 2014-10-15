#!/bin/bash 
# --- Check AliEn token ----------------------------------------------
uid=`id -u`
genv_file=/tmp/gclient_env_${uid}

if test ! -f ${genv_file} ; then 
    echo "No such file: ${genv_file}, please do alien-token-init" \
        >/dev/stderr
    exit 1
fi
. ${genv_file}
alien-token-info | grep -q "Token is still valid"
if test $? -ne 0 ; then 
    echo "Token not valid, please re-new" > /dev/stderr 
    exit 1
fi

run=118506
nev=1
if test x$1 != x ; then run=$1 ; fi 
files="AOD.C	\
	AODConfig.C	\
	Check.C		\
	Config.C	\
	GRP.C		\
	QA.C		\
	QAConfig.C	\
	Reconstruct.C	\
	Simulate.C	\
	Tag.C		\
	simrun.sh	\
	$ALICE_ROOT/OADB/PWGLF/FORWARD/CORRECTIONS/data/fmd_corrections.root"

rm -rf test
mkdir -p test
for i in $files ; do 
    cp -v $i test/`basename $i` 
done 

(cd test && ../run.sh --run $run  --event $nev --qa --aod $@)
