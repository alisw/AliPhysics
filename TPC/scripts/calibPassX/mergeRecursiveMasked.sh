#
# recursive merging
# 
maxMerge=$1
queue="$2 -m batch_dgrid2 -E /u/miranov/preexec.sh"
mask=$3
output=$4
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
     echo    bsub -q $queue  -oo outMerge.log $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"$output\",\"$mask\",\"riend\"\);
     cat calib.list
     bsub -q $queue  -oo outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"$output\",\"$mask\",\"riend\"\)
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
#echo bsub -q $queue  -oo outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"$output\","\$mask\",\"riend\"\)
bsub -q $queue  -oo outMerge.log aliroot -b -q  $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"$output\",\"$mask\",\"riend\"\)

