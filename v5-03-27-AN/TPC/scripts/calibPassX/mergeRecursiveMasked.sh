#
# recursive merging
# 
maxMerge=$1
queue="$2"
mask=$3
output=$4
reject=$5
#
counter=0;
counter2=0;
wdir=`pwd`
rm -rf merge*
mkdir merge$counter2 
cd merge$counter
for a in `cat ../calib.list`; do
   let counter=counter+1;
   echo $counter $counter2
   echo $a >>calib.list
   if [ $counter -gt $maxMerge ] ; then
     echo    bsub -q $queue  -oo outMerge.log $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"$output\",\"$mask\",\"$5\"\);
     cat calib.list
     bsub -q $queue  -oo outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"$output\",\"$mask\",\"$5\"\)
     let counter2=counter2+1;
     let counter=0;
     cd $wdir
     mkdir  merge$counter2 
     cd merge$counter2
     if [ -e calib.list ]; then
        rm calib.list
     fi;
   fi;  
done;

bsub -q $queue  -oo outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"$output\",\"$mask\",\"$5\"\)

