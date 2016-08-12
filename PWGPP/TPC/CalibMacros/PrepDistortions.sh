#!/bin/bash

cdbPath="local://"

if [[ $# < 2 ]]; then
 echo "Append distortion Cheb.param to AliTPCDcalibRes with already processed correction"
 echo "and extract distortion CDB object to $cdbPath"
 echo ""
 echo "Usage: ./${0##*/} <path_to_timebins_topdir> <run> [produce_test_tree]"
 echo "1) path_to_timebins_topdir is directory containing the <tmin>_<tmax>_<run> subdirectories"
 echo "   If this directory on the grid, precede it by alien://"
 echo "2) run number"
 echo "3) optional non-empty argument to request extraction of control tree with corrections and distortions" 
 echo "Example: ./${0##*/} alien:///alice/data/2015/LHC15n/000244540/cpass0_pass1/ResidualMerge/TPCSPCalibration 244540"
 exit 1
fi


clbname="alitpcdcalibres.root"

macroInv="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/PrepDistortion.C"

cdbUpdScript="$ALICE_PHYSICS/PWGPP/TPC/CalibMacros/createOCDBTPCSPDistCalib.sh"

pathTB=$1
run=$2

[[ -z $3 ]] && testTree=0 || testTree=1

echo "data path  $pathTB"
echo "run:       $run"

if [[ $pathTB == 'alien://'* ]] ; then
    gridSrc='1'
    pathTB=${pathTB/"alien://"/""}
    lst=( `alien_find "$pathTB" "*_*_000*/$clbname" | grep "$pathTB" | perl -p -e 's/$pathTB/alien:\/\/$pathTB/'` )
else 
    lst=( `find $pathTB -path "*_*_000*/$clbname"` )
fi


nchunks="${#lst[@]}"

echo "Found $nchunks chunks to process:" 
echo "${lst[@]}"

curdir=`pwd`

outList="resList.txt"
[[ -e $outList ]] && rm $outList 

locMacroInv=$(basename "$macroInv")
[[ ! -f "$locMacroInv" ]] && cp $macroInv ./ 
[[ ! -f "$locMacroInv" ]] && echo "did not find $locMacroInv" && exit 1


for chunk in ${lst[@]} 
do  
    echo "Processing $chunk"
    # extract the directory name of the timebin
    dir="$(dirname $chunk)"
    dirTB="$(basename $dir)"
    [[ -n $dirTB ]] && mkdir $dirTB ; cd $dirTB
    echo aliroot -b -q "$curdir/$locMacroInv"\(\"$chunk\",$testTree\)
    aliroot -b -q "$curdir/$locMacroInv"\(\"$chunk\",$testTree\)
    [[ -e $clbname ]] && echo `pwd`/"$clbname" >> "$curdir/$outList" || echo "Processing failed in `pwd`"
    cd $curdir
done
 
# create distortion ocdb object
corr="0"
dist="1"
$cdbUpdScript inputFileList=$outList startRun=$run endRun=$run ocdbStorage=$cdbPath corr=$corr dist=$dist
