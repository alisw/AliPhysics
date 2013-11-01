            
source $1 
aliroot -b -q $ALICE_ROOT/TPC/Upgrade/macros/fastToyMCMI.C+
exit;


# Example usage local 
# jobs to be submitted form the lxb1001 or lxb1002
#(here we have 80 nodes and user disk)
#
export baliceTPC=/u/miranov/.baliceTPC
export flucPath=$HOME/AliRoot/TPCdev/TPC/Upgrade/macros/
export batchCommand="qsub -cwd  -V "

for idir in {0..40}; do  
   $batchCommand    -o  drawFlucBin.log  $flucPath/fastToyMCMI.sh  $baliceTPC  7  100000 10000 0 0
done;

wdir=`pwd`
for a in `ls -d dir*`; do
    cd $wdir/$a
    rm localBins.root
    $batchCommand    -o  drawFlucBin.log  $flucPath/spaceChargeFluctuation.sh $baliceTPC  7  100000 10000 0 0
    cd $wdir
done;

