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

EDIR=TPCcalib

mkdir $EDIR


SRC=$ALICE_ROOT/TPC/TPCcalib
echo Source $SRC  
echo EDIR $EDIR

cp $SRC/Ali*.h                $EDIR
cp $SRC/Ali*.cxx              $EDIR
cp $SRC/TPCcalibLinkDef.h     $EDIR
cp $SRC/Makefile*             $EDIR
cp $SRC/libTPCcalib.pkg2      $EDIR  


mkdir $EDIR/PROOF-INF
cd $EDIR/PROOF-INF


cp $SRC/BUILD.sh .
cp $SRC/SETUP.C  .


chmod 755 BUILD.sh

cd ../..

tar zcvf TPCcalib.par $EDIR


exit 0
