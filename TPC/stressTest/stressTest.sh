# Run stress test on the batch system
# All run*.sh scripts in the $ALICE_ROOT test macro invoked
#
# Parameters:
# 1         - output prefix
# 2         - submit command
# 
# Run example:
# $ALICE_ROOT/test/stressTest/stressTest.sh /d/alice12/miranov/streeTest/ "bsub -q proof"
# 

outdir=$1/$ALICE_LEVEL/$ALICE_TARGET
submitcommand=$2
echo _____________________________________________________________
echo _____________________________________________________________
echo _____________________________________________________________
echo
echo outdir        $outdir
echo subitcommand  $submitcommand
mkdirhier $outdir
ls -al    $outdir
echo
echo _____________________________________________________________
echo _____________________________________________________________
echo _____________________________________________________________

#
# Loop over all run*sh macros
#
for tmacro in `ls $ALICE_ROOT/test/*/run*.sh` ; do
#
dname=`dirname $tmacro`
sname=`basename $dname`
workdir=$outdir/$sname
echo $sname $tmacro;
mkdirhier $workdir
cp $dname/* $workdir/
cd $workdir
rm *.root
echo $submitcommand -oo $workdir/out.log --eo $workdir/err.log $tmacro
$submitcommand -oo $workdir/out.log -eo $workdir/err.log $tmacro
cd $outdir;
done;





