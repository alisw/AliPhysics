#  
# shell script to submit small jobs for the current makeCurrentCorrection.C root macro
# used for the TPC TDR studies 
#         
source $1 
aliroot -b -q $flucPath/NimStyle.C $ALICE_ROOT/TPC/Upgrade/macros/makeCurrentCorrection.C+\($2,$3,$4\)
exit;

#
# Example usage local jobs to be submitted form the lxb1001 or lxb1002
#(here we have 80 nodes and user disk)
#

makeEnvLocal(){
#
#
#
    export baliceTPC=/u/miranov/.baliceTPC
    export flucPath=$HOME/AliRoot/git/TPCdev/TPC/Upgrade/macros/
    export batchCommand="qsub -cwd  -V "
}

submitCurrentExtractionJobs(){
#
# Action
# 1.) submit current extraction jobs
#
    wdir=`pwd`
    counter=0;
    for a in `ls -d dir*`; do
	cd $wdir/$a
	rm localCurent.root current.log 
	$batchCommand    -o  current.log  $flucPath/makeCurrentCorrection.sh $baliceTPC  0  20000 $counter
	let counter=counter+1
	cd $wdir
    done;
}

submitSmoothingJobs(){
#
# Action:
# 2.) Submit smoothing jobs
#
    for imap  in {0..40}; do
	mkdir /hera/alice/miranov/SpaceCharge/Smoothing2D_2/test$imap
	cd  /hera/alice/miranov/SpaceCharge/Smoothing2D_2/test$imap
	ln -sf /hera/alice/wiechula/Upgrade/LUTs_fluctuation_eps20/MeasuredResidual2D/residualMap.$imap.root MeasureResidual.root
	ln -sf /hera/alice/wiechula/Upgrade/LUTs_fluctuation_eps20/RealResidualScaled/residualMap.$imap.root RealResidualScaled.root
	cd  /hera/alice/miranov/SpaceCharge/Smoothing2D_2/test$imap
	$batchCommand    -o  smoothing.log  $flucPath/makeCurrentCorrection.sh $baliceTPC  1 2000 80 $imap
    done;
}

