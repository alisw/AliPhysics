#
# run PassX recosntruction/calibration jobs
# 
# Needed in order to run the jobs locally
#
#    1  - raw data file name
#    2  - number of events to be processed
#    3  - run number 


echo Input file $1
echo Run number $3
echo Events to reconstruct/calibrate $2

outdir=`pwd`
tmpdir=/tmp/$USER/`echo $outdir | sed s_/_x_g` 
mkdir $tmpdir
#
cp  $outdir/*.C   $tmpdir/
cp  $outdir/*.sh  $tmpdir/
cp  $outdir/AliESDs.root        $tmpdir/
cp  $outdir/AliESDfriends.root  $tmpdir/

cd $tmpdir
#
#
#
echo tmpdir = `pwd`
echo outdir = $outdir
ls

echo runPassX.sh   $1  $2  $3
runPassX.sh        $1  $2  $3

ls -alrt 
cp -r $tmpdir/*      $outdir/
rm -rf $tmpdir


