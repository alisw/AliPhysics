#!/bin/sh

# raw2sdigits.sh rawfile [esdfile]
#
# convert a (real) raw data file into SDigits suitable for embedding.
# if esdfile is present (should correspond to the raw data file ;-) )
# then zvertex is taken from it

rawfile=$1
esdfile=""

if [ $# -ge 2 ]; then
 esdfile=$2
fi

echo "rawfile=$rawfile esdfile=$esdfile"

aliroot -l -b << EOF

    AliSimulation sim("$ALICE_ROOT/MUON/Config.C");
    gAlice->SetConfigFunction("Config(\"\", \"param\", \"AliMUONDigitStoreV2R\",kTRUE);");
    AliCDBManager *cdbm = AliCDBManager::Instance();
    AliRawReader* reader = AliRawReader::Create(gSystem->ExpandPathName("$rawfile"));
    reader->NextEvent();
    Int_t runNumber = reader->GetRunNumber();
    delete reader;
    cout << "SETTING RUN NUMBER TO " << runNumber << endl;
    cdbm->SetRun(runNumber);
    cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
    cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
    sim.SetMakeSDigits("MUON");  

    sim.ConvertRaw2SDigits("$rawfile","$esdfile");

EOF
