#!/bin/bash

mc=0
rew=0
class=MakeTrackletTrain
run=245064
daper=LHC15o
mcper=LHC15k1_plus21
aliph="last,normal"
aliph="vAN-20160207-1"
dadir="data"
mcdir="sim";
dapat="pass_lowint_firstphys/*/AliESDs.root"
mcpat="*/AliESDs.root"
daopt="&split=50"
mcopt="&split=100&mc"

case $1 in
    --*|-*) ;;
    mc) mc=1 ; shift ;;
    rmc) mc=1 ; rew=1 ; shift ;; 
    *) ;;
esac

if test $mc -gt 0 ; then
    per=$mcper
    dir=$mcdir
    pth=$mcpth
    pat=$mcpat
    opt=$mcopt
else 
    per=$daper
    dir=$dadir
    pth=$dapth
    pat=$dapat
    opt=$daopt
fi
if test $rew -gt 0 ; then
    more="--reweight=pt,pid,str"
    xtra="_reweighted"
fi
yer=`echo $per | sed 's/LHC\(..\).*/20\1/'` 
pth="alien:///alice/${dir}/${yer}/${per}"
url="${pth}/?run=${run}&pattern=${pat}&aliphysics=${aliph}${opt}#esdTree"
nme="${per}_${run}_tracklets${xtra}"

opts=(--class=$class \
	     --name="$nme" \
	     --cent="V0M" \
	     --cent-bins="0-5-10-20-30-40-50-60-70-80" \
	     --cent-oadb="" \
	     --create-inj \
	     --events="-1" \
	     --trig="V0AND" \
	     --url="${url}" \
	     --verbose="0" \
	     $more)

# echo "Will do runTrain ${opts[@]} $@"
# exit

rm -rf $nme
auser=`alien_whoami | sed -e 's/^ *//' -e 's/ *$//'` 
adir="/alice/cern.ch/user/c/${auser}"
alien_rmdir ${adir}/${nme}

echo "Will do runTrain ${opts[@]} $@"
runTrain "${opts[@]}" $@

#
# EOF
#

