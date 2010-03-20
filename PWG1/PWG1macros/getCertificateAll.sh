# Pseudo code to test the all components in the train and make report


#
# 1. copy all task ans shell scirpts
#
cp /u/miranov/AliRoot/trunk/PWG1/PWG1macros/*.* .
 

#
# 2. Submit jobs for each macro - in separate directory
#
bqueue=alice-t3_8h
workdir=`pwd`
for amacro in `ls  Add*.C`; do
    mkdirhier $workdir/test$amacro
    cd $workdir/test$amacro
    cp $workdir/$amacro .
    cp $workdir/runPWG1Train.C .
    cp $workdir/*.sh  .
    cp $workdir/esd.list .
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
