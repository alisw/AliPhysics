if [ "0$1" == "0" ]; then
    MYFILE=raw.root
else
    MYFILE=$1
fi

MYOCDB=local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2018/OCDB/

rm -f galice.root
rm -f *.QA.*.root
rm -Rf testDir

MYPREPENDFILE=test.C

MYFIRSTEVENT=-1
MYLASTEVENT=-1

stdbuf -o0 aliroot -l -q -b HLTVZeroCalibTest.C $ALICE_ROOT/HLT/exa/recraw-local.C"(\"$MYFILE\",\"$MYOCDB\", $MYFIRSTEVENT, $MYLASTEVENT, \"HLT\", \"chains=RootWriter ignore-hltout loglevel=0x7c\")" | tee runhltvzero.log
