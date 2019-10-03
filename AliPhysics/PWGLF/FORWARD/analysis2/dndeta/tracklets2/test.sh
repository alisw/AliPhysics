#!/bin/bash

type=dt
case x$1 in
    x--*) ;;
    x-*) ;;
    x) ;;
    x*) type=$1 ; shift ;;
esac

run=244918
class="TrackletdNdetaTrain"
name="LHC15o_${run}_midRapidity_data"
idx=/data/alice/data/pbpb/LHC15o/000244918_fp/index.root
opts=
if test "x$type" = "xmc" ; then
    idx=/data/alice/data/pbpb/LHC15o/000244918_mc/index.root
    opts="&mc"
    name="LHC15o_${run}_midRapidity_sim"
fi 
url="lite:///${idx}?workers=10${opts}#esdTree"

opts=(--class="$class" \
	     --name="${name}" \
	     --date="none" \
	     --cent="V0M" \
	     --cent-bins="0-5-10-20-30-40-50-60-70-80" \
	     --eta-bins="r16:2" \
	     --ipz-bins="-15:-12.5:-10:-8:-6:-4:-2:0:2:4:6:8:10:12.5:15" \
	     --ipz-bins="u:15" \
	     --reconstruct="nor,inj" \
	     --ocdb \
	     --trig="V0AND" \
	     --url="${url}" \
	     --verbose="0" \
	     --no-link \
	     ${more})

rm -rf $name
echo "runTrain ${opts[@]} $@"
runTrain ${opts[@]} $@
