#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "TCollection.h"
#include "TMap.h"
#include "TKey.h"
#include "TString.h"
#include "TObjString.h"
#include "TFile.h"
#include "TSystem.h"

void MakeSnapshot(Int_t run, const char* defStorage, TMap* specStorages, const char* snapshotFileName)
{
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage(defStorage);
    cdb->SetRun(run);

    TIter iter(specStorages->GetTable());
    TPair *pair = 0;
    while((pair = dynamic_cast<TPair*> (iter.Next()))){
	TObjString* caltype = dynamic_cast<TObjString*> (pair->Key());
	TObjString* specstor= dynamic_cast<TObjString*> (pair->Value());
	if (caltype && specstor)
	    //TString calType = caltype->GetString();
	    //TString specStor = specstor->GetString();
	    //cdb->SetSpecificStorage(calType.Data(),specStor.Data());
	    cdb->SetSpecificStorage(caltype->GetString().Data(),specstor->GetString().Data());
	else
	    //AliFatal("Error reading info for specific storage")
	    Printf("Error reading info for specific storage");
    }

    // ********************************** GRP ******************************************
    cdb->Get("GRP/CTP/Config");
    cdb->Get("GRP/Calib/LHCClockPhase");
    cdb->Get("GRP/GRP/Data");
    cdb->Get("GRP/Align/Data");
    cdb->Get("GRP/Calib/MeanVertexSPD");
    cdb->Get("GRP/Calib/MeanVertex");
    cdb->Get("GRP/Calib/MeanVertexTPC");
    cdb->Get("GRP/Calib/CosmicTriggers");
    cdb->Get("GRP/CTP/Scalers");
    cdb->Get("GRP/CTP/CTPtiming");
    cdb->Get("GRP/CTP/TimeAlign");
    cdb->Get("GRP/GRP/LHCData");
    cdb->Get("GRP/Calib/RecoParam");

    // ********************************** ALL ******************************************
    TString detStr = ("ITS TPC TRD TOF PHOS HMPID EMCAL MUON ZDC PMD T0 VZERO");
    //TString detStr = ("ITS MUON TPC");
    TObjArray *arr = detStr.Tokenize(' ');
    for (Int_t iDet=0; iDet<arr->GetEntries(); iDet++) {
	TObjString *detOStr = dynamic_cast<TObjString*>(arr->At(iDet));
	AliCDBManager::Instance()->GetAll(Form("%s/Calib/*",detOStr->GetString().Data()));
	AliCDBManager::Instance()->Get(Form("%s/Align/Data",detOStr->GetString().Data()));
    }

    // ******************************** TRIGGER ****************************************
    // Temporary fix - one has to define the correct policy in order
    // to load the trigger OCDB entries only for the detectors that
    // in the trigger or that are needed in order to put correct
    // information in ESD
    AliCDBManager::Instance()->GetAll("TRIGGER/*/*");

    // ********************************** HLT ******************************************
    // cdb->Get("HLT/ConfigHLT/esdLayout");
    // cdb->Get("HLT/Calib/StreamerInfo");

    TMap* entriesMap = const_cast<TMap*>(cdb->GetEntryCache());
    Printf("\nentriesMap has %d entries!\n", entriesMap->GetEntries());

    TList* entriesList = const_cast<TList*>(cdb->GetRetrievedIds());
    Printf("\nentriesList has %d entries!\n", entriesList->GetEntries());

    //TString filename(TString::Format("CDBsnapshot_Run%d.root",run));
    TString filename(snapshotFileName);
    TFile *f = new TFile(filename.Data(),"recreate");
    f->cd();
    f->WriteObject(entriesMap,"entriesMap");
    f->WriteObject(entriesList,"entriesList");
    f->Close();
    entriesMap->SetOwnerKeyValue(kFALSE,kFALSE);
    entriesList->SetOwner(kFALSE);
}

void MakeCDBSnapshot()
{
    //TMap *specMap = new TMap(10);
    TMap *specMap = new TMap(30);
    specMap->SetName("mapOfSpecificStorages");
    specMap->SetOwner(1);
    specMap->Add(new TObjString("ITS/Align/Data"), new TObjString("alien://folder=/alice/simulation/2008/v4-15-Release/Residual"));

    specMap->Add(new TObjString("MUON/Align/Data"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("MUON/Calib/Gains"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    //      MTR
    specMap->Add(new TObjString("MUON/Calib/GlobalTriggerCrateConfig"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Full"));
    specMap->Add(new TObjString("MUON/Calib/LocalTriggerBoardMasks"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Full"));
    specMap->Add(new TObjString("MUON/Calib/RegionalTriggerConfig"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Full"));
    specMap->Add(new TObjString("MUON/Calib/TriggerEfficiency"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Full"));
    // TPC (23 total) 
    specMap->Add(new TObjString("TPC/Calib/PadGainFactor"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/TimeGain"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/GainFactorDedx"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/PadTime0"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Distortion"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/PadNoise"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/PadNoise"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Pedestals"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Temperature"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Parameters"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/ClusterParam"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/AltroConfig"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Pulser"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Pulser"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/CE"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Mapping"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Correction"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Align/Data"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Goofie"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/TimeDrift"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/Raw"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/QA"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/"));
    specMap->Add(new TObjString("TPC/Calib/HighVoltage"), new TObjString("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/")); 


    MakeSnapshot(139517,"alien://folder=/alice/data/2010/OCDB",specMap,"cdbSnapshot.root");

}
