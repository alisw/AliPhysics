// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TTimeStamp.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include "AliTimeStamp.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliADCalibData.h"


void AnalyzeCTPRecords(Int_t rn, Int_t year, const TString *classNames, Int_t nClasses) {
  Printf("nClasses=%d", nClasses);;

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(Form("alien://folder=/alice/data/%d/OCDB", year));

  TTree *TS = new TTree;
  TS->SetName("TS");

  TMatrixD ns(nClasses, 6);
  TMatrixD offset(nClasses,6);
  TTimeStamp *timeStamp = NULL;
  TS->Branch("TimeStamp", "TTimeStamp", &timeStamp, 32000, 1);
  for (Int_t i=0; i<nClasses; ++i) {
    TObjArray *os = classNames[i].Tokenize("-");
    TS->Branch(Form("%s_%s", os->At(0)->GetName(), os->At(1)->GetName()), ns.GetMatrixArray()+6*i, "L0b/D:L0a:L1b:L1a:L2b:L2a");
  }
  man->SetRun(rn);

  AliCDBEntry *entry = man->Get("GRP/CTP/Config");
  AliTriggerConfiguration *triggerConfig  = (AliTriggerConfiguration*)entry->GetObject();

  entry = man->Get("GRP/CTP/Scalers");
  AliTriggerRunScalers    *triggerScalers = (AliTriggerRunScalers*)entry->GetObject();

  TArrayI classIndex(nClasses);
  for (Int_t i=0; i<nClasses; ++i) {
    classIndex[i] = triggerConfig->GetClassIndexFromName(classNames[i]);
    if (classIndex[i] < 0) {
      Printf("class '%s' not found", classNames[i].Data());
      return;
    }
  }

  const TObjArray *a = triggerScalers->GetScalersRecords();
  for (Int_t j=0, m=a->GetEntries(); j<m; ++j) {
    AliTriggerScalersRecord *scalersRecord = (AliTriggerScalersRecord*)a->At(j);
    const AliTimeStamp *tt = scalersRecord->GetTimeStamp();
    *timeStamp = TTimeStamp(time_t(tt->GetSeconds()), 1000*tt->GetMicroSecs());

    Printf("time: %s (%d %d) ", timeStamp->AsString(), tt->GetSeconds(), tt->GetMicroSecs());

    for (Int_t i=0; i<nClasses; ++i) {
      const AliTriggerScalers *scalers = scalersRecord->GetTriggerScalersForClass(classIndex[i]);
      if (j == 0) {
	offset(i,0) = scalers->GetLOCB();
	offset(i,1) = scalers->GetLOCA();
	offset(i,2) = scalers->GetL1CB();
	offset(i,3) = scalers->GetL1CA();
	offset(i,4) = scalers->GetL2CB();
	offset(i,5) = scalers->GetL2CA();
      }
      ns(i,0) = scalers->GetLOCB() - offset(i,0);
      ns(i,1) = scalers->GetLOCA() - offset(i,1);
      ns(i,2) = scalers->GetL1CB() - offset(i,2);
      ns(i,3) = scalers->GetL1CA() - offset(i,3);
      ns(i,4) = scalers->GetL2CB() - offset(i,4);
      ns(i,5) = scalers->GetL2CA() - offset(i,5);
      Printf("%s %.0f", classNames[i].Data(), ns(i,0));
    }
    TS->Fill();
  }

  TFile::Open(Form("root/scalers_%d.root", rn), "RECREATE");
  TS->Write();
  gFile->Write();
  gFile->Close();
}
