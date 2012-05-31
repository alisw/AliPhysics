#! /bin/sh
#
# Make the event.par file to be used to analyse Event objects with PROOF.
#
# Usage: sh make_event_par.sh
#
# Creates the PAR file "event.par" which can be used in PROOF via the
# package manager like:
#   gProof->UploadPackage("event.par")
#   gProof->EnablePackage("event")
#
# Command to check that package is active and that libEvent.so is loaded:
#   gProof->ShowPackages()
#   gProof->ShowEnabledPackages()
#   gProof->Exec("gSystem->ListLibraries()")
#

EDIR=TPCcalibPar

mkdir $EDIR


SRC=$ALICE_ROOT/TPC
echo Source $SRC  
echo EDIR $EDIR
cp $SRC/AliTPCFitPad*               $EDIR
cp $SRC/AliTPCCal*.h                $EDIR
cp $SRC/AliTPCCal*.cxx              $EDIR
cp $SRC/AliTPCcal*.h                $EDIR
cp $SRC/AliTPCcal*.cxx              $EDIR
cp $SRC/AliTPCSel*.cxx              $EDIR
cp $SRC/AliTPCSel*.h                $EDIR
cp $SRC/AliTPC*Ana*.*               $EDIR
cp $SRC/TPCcalibLinkDef.h           $EDIR
cp $SRC/Makefile.Calib              $EDIR/Makefile
cp $SRC/Makefile.arch.Calib         $EDIR/Makefile.arch
cp $SRC/libTPCcalib.pkg             $EDIR  


mkdir $EDIR/PROOF-INF
cd $EDIR/PROOF-INF


cp $SRC/BUILDcalib.sh BUILD.sh
cp $SRC/SETUPcalib.C  SETUP.C 


chmod 755 BUILD.sh

cd ../..

tar zcvf TPCcalibPar.par $EDIR


exit 0
