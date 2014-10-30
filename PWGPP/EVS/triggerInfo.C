#ifndef __CINT__
#include "TMath.h"
#include "AliCDBManager.h"
#include "AliTriggerScalers.h"
#include "AliTriggerRunScalers.h"
#include "AliTimeStamp.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerConfiguration.h"
#include "AliLHCData.h"
#include "AliTriggerClass.h"
#include "AliTriggerBCMask.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#endif

UInt_t dif(UInt_t stop, UInt_t start){
  UInt_t d;
  if(stop >= start) d=stop-start;
  else d = stop+(0xffffffff-start)+1;
  return d;
};

Int_t triggerInfo(Int_t run, TString refClassName, Double_t* par){
  AliCDBManager* man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("local:///data/alice/OCDB");
  man->SetRun(run);
  AliLHCData* lhc = (AliLHCData*) man->Get("GRP/GRP/LHCData")->GetObject();
  par[0] = lhc->GetFillNumber();
  
  AliTriggerConfiguration* cfg = (AliTriggerConfiguration*) man->Get("GRP/CTP/Config")->GetObject();

  // Get scalers 
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*) man->Get("GRP/CTP/Scalers")->GetObject();
  Int_t nEntries = scalers->GetScalersRecords()->GetEntriesFast();
  
  // Get SOR and EOR scaler records
  AliTriggerScalersRecord* record1 = scalers->GetScalersRecord(0);
  AliTriggerScalersRecord* record2 = scalers->GetScalersRecord(nEntries-1);
  if (!record1) { printf("Null pointer to start scalers record\n"); return 1; }
  if (!record2) { printf("Null pointer to  stop scalers record\n"); return 1; }

  // Extract SOR and EOR times
  const AliTimeStamp* stemp1 = record1->GetTimeStamp();
  const AliTimeStamp* stemp2 = record2->GetTimeStamp();
  if (!stemp1) { printf("Null pointer to start timestemp\n"); return 1; }
  if (!stemp2) { printf("Null pointer to stop timestemp\n");  return 1; }
  Double_t duration = stemp2->GetSeconds()-stemp1->GetSeconds();
  par[1] = duration;
  if (TMath::Abs(duration)<1) return 2;

  // Extract SOR and EOR trigger counts
  Int_t classid = cfg->GetClassIndexFromName(refClassName);
  const AliTriggerScalers* scaler1 = record1->GetTriggerScalersForClass(classid);
  const AliTriggerScalers* scaler2 = record2->GetTriggerScalersForClass(classid);
  if (!scaler1) { printf("Null pointer to start scalers for reference class\n"); return 1; }
  if (!scaler2) { printf("Null pointer to stop scalers for reference class\n");  return 1; }
  UInt_t l0b = dif(scaler2->GetLOCB(),scaler1->GetLOCB());
  par[2] = l0b;

  // Get number of colliding bunches per orbit
  Double_t orbitRate = 11245.; // Hz
  Double_t nBCsPerOrbit = -1;
  if (refClassName.Contains("1B-ABCE-")){
    nBCsPerOrbit = lhc->GetNInteractingBunchesMeasured();
    Printf("Number of BCs from LHC data=%i",nBCsPerOrbit);
    if (nBCsPerOrbit<0) {
      Int_t emptyclassid = cfg->GetClassIndexFromName("CBEAMB-ABCE-NOPF-ALL");
      if (emptyclassid<0) return 3;
      const AliTriggerScalers* emptyScaler1 = record1->GetTriggerScalersForClass(emptyclassid);
      const AliTriggerScalers* emptyScaler2 = record2->GetTriggerScalersForClass(emptyclassid);
      if (!scaler1) { printf("Null pointer to start scalers for reference class\n"); return 1; }
      if (!scaler2) { printf("Null pointer to stop scalers for reference class\n");  return 1; }
      UInt_t l0bempty = dif(emptyScaler2->GetLOCB(),emptyScaler1->GetLOCB());
      if (l0bempty==0) return 4;
      nBCsPerOrbit = l0bempty/orbitRate/duration;
    }
  }
  else {
    // Extract number of bunches per orbit
    AliTriggerClass* cl = cfg->GetTriggerClass(classid);
    AliTriggerBCMask* mask = cl->GetBCMask();
    nBCsPerOrbit = mask->GetNUnmaskedBCs();
  }
  par[3] = nBCsPerOrbit;
  
  Double_t totalBCs = duration*orbitRate*nBCsPerOrbit;
  par[4] = -TMath::Log(1-l0b/totalBCs); // mu
  return 0;
}
