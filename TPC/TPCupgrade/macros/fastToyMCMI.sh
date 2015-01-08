            
source $1 
aliroot -b -q $ALICE_ROOT/TPC/Upgrade/macros/fastToyMCMI.C+\($2,$3,$4)
exit;


# Example usage local 
# jobs to be submitted form the lxb1001 or lxb1002
#(here we have 80 nodes and user disk)
#
export baliceTPC=/u/miranov/.baliceTPC
export flucPath=$HOME/AliRoot/TPCdev/TPC/Upgrade/macros/
export batchCommand="qsub -cwd  -V "

wdir=`pwd`
for idir in {0..20}; do  
    mkdir $wdir/dir$idir
    cd  $wdir/dir$idir
    ln -sf ../geometry.root .
    rm testdEdxResolution.root
   $batchCommand    -o  fastToyMCMI.log  $flucPath/fastToyMCMI.sh  $baliceTPC  1 200 0 $idir
   
done;

