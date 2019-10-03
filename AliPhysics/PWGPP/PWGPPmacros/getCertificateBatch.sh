#
# Pseudo code to test the all components in the train and make report on batch farm
#

#
# 1. copy all task ans shell scirpts
#
cp /u/miranov/AliRoot/trunk/PWGPP/PWGPPmacros/*.* .
 

#
# 2. Submit jobs for each macro - in separate directory
#
bqueue=alice-t3_8h
workdir=`pwd`
rm -rf test*
for fmacro in `cat ConfigTask.txt`; do
    amacro=`basename $fmacro`
    dname=`echo $workdir/test$amacro | sed s_.C__`
    mkdirhier $dname
    cd $dname
    #
    cp $workdir/runPWGPPTrain.C .
    cp $workdir/*.sh  .
    cp $workdir/esd.list .
    cp $workdir/ConfigTask.txt .
    echo bsub -q $bqueue getCertificate.sh $amacro
    bsub -q $bqueue getCertificate.sh $amacro esd.list
    cd $workdir
done

#
# 3. Wait
#
echo name/C:time/C > summaryTime.txt
for a in `ls */summary.log` ; do echo $a `cat $a | grep SysInfoTime` >>summaryTime.txt ; done

echo name/C: time/C > sumaryMem.txt
for a in `ls */summary.log` ; do echo $a `cat $a | grep SysInfoMem` >>summaryMem.txt ; done
