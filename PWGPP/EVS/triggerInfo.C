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

TObjArray GetClasses(Int_t run, TString ocdbStorage, ULong64_t* L0B, ULong64_t* L0A, ULong64_t* L1B, ULong64_t* L1A, ULong64_t* L2B, ULong64_t* L2A){
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbStorage.Data());
  man->SetRun(run);
  // Get scalers 
  AliTriggerConfiguration* cfg = (AliTriggerConfiguration*) man->Get("GRP/CTP/Config")->GetObject();
  if (!cfg) { printf("No GRP/CTP/Config object for run %i\n",run); return TObjArray(); }

  TObjArray classes = cfg->GetClasses();
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*) man->Get("GRP/CTP/Scalers")->GetObject();
  if (!scalers) { printf("No GRP/CTP/Scalers object for run %i\n",run); return TObjArray(); }
  Int_t nEntries = scalers->GetScalersRecords()->GetEntriesFast();
  
  for (Int_t r=0;r<nEntries-1;r++){
    // Get SOR and EOR scaler records
    AliTriggerScalersRecord* record1 = scalers->GetScalersRecord(r);
    AliTriggerScalersRecord* record2 = scalers->GetScalersRecord(r+1);
    if (!record1) { printf("Null pointer to scalers record\n"); return TObjArray(); }
    if (!record2) { printf("Null pointer to scalers record\n"); return TObjArray(); }
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      // Extract SOR and EOR trigger counts
      Int_t classId = cfg->GetClassIndexFromName(classes.At(i)->GetName());
      const AliTriggerScalers* scaler1 = record1->GetTriggerScalersForClass(classId);
      const AliTriggerScalers* scaler2 = record2->GetTriggerScalersForClass(classId);
      if (!scaler1) { printf("Null pointer to scalers for class\n"); return TObjArray(); }
      if (!scaler2) { printf("Null pointer to scalers for class\n"); return TObjArray(); }
      L0B[i] += dif(scaler2->GetLOCB(),scaler1->GetLOCB());
      L0A[i] += dif(scaler2->GetLOCA(),scaler1->GetLOCA());
      L1B[i] += dif(scaler2->GetL1CB(),scaler1->GetL1CB());
      L1A[i] += dif(scaler2->GetL1CA(),scaler1->GetL1CA());
      L2B[i] += dif(scaler2->GetL2CB(),scaler1->GetL2CB());
      L2A[i] += dif(scaler2->GetL2CA(),scaler1->GetL2CA());
    }
  }
  return classes;
}


Int_t triggerInfo(Int_t run, TString refClassName, TString ocdbStorage, TString &activeDetectors, Double_t* par){
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbStorage.Data());
  man->SetRun(run);
  if (run!=189694) { // No GRP/GRP/LHCData for this run
    AliLHCData* lhc = (AliLHCData*) man->Get("GRP/GRP/LHCData")->GetObject();
    if (!lhc) { printf("No GRP/GRP/LHCData object for run %i\n",run); return 1; }
    par[0] = lhc->GetFillNumber();
  }
  
  AliTriggerConfiguration* cfg = (AliTriggerConfiguration*) man->Get("GRP/CTP/Config")->GetObject();
  if (!cfg) { printf("No GRP/CTP/Config object for run %i\n",run); return 1; }
  activeDetectors = cfg->GetActiveDetectors().Data();

  // Get scalers 
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*) man->Get("GRP/CTP/Scalers")->GetObject();
  if (!scalers) { printf("No GRP/CTP/Scalers object for run %i\n",run); return 1; }
  Int_t nEntries = scalers->GetScalersRecords()->GetEntriesFast();

  Double_t run_duration   = 0;
  ULong64_t l0b      = 0;
  ULong64_t l0bempty = 0;
  Int_t classId      = cfg->GetClassIndexFromName(refClassName);
  Int_t emptyclassid = cfg->GetClassIndexFromName("CBEAMB-ABCE-NOPF-ALL");
  for (Int_t r=0;r<nEntries-1;r++){
    // Get SOR and EOR scaler records
    AliTriggerScalersRecord* record1 = scalers->GetScalersRecord(r);
    AliTriggerScalersRecord* record2 = scalers->GetScalersRecord(r+1);
    if (!record1) { printf("Null pointer to scalers record\n"); return 2; }
    if (!record2) { printf("Null pointer to scalers record\n"); return 2; }
    const AliTimeStamp*      stamp1  = record1->GetTimeStamp();
    const AliTimeStamp*      stamp2  = record2->GetTimeStamp();
    const AliTriggerScalers* scaler1 = record1->GetTriggerScalersForClass(classId);
    const AliTriggerScalers* scaler2 = record2->GetTriggerScalersForClass(classId);
    if (!stamp1 ) { printf("Null pointer to start timestamp\n");   return 2; }
    if (!stamp2 ) { printf("Null pointer to stop timestamp\n");    return 2; }
    if (!scaler1) { printf("Null pointer to scalers for class %s\n",refClassName.Data()); return 2; }
    if (!scaler2) { printf("Null pointer to scalers for class %s\n",refClassName.Data()); return 2; }
//    run_duration += dif(stamp2->GetSeconds(),stamp1->GetSeconds());
    l0b          += dif(scaler2->GetLOCB()  ,scaler1->GetLOCB());
    if (emptyclassid<0) continue;
    const AliTriggerScalers* emptyScaler1 = record1->GetTriggerScalersForClass(emptyclassid);
    const AliTriggerScalers* emptyScaler2 = record2->GetTriggerScalersForClass(emptyclassid);
    if (!emptyScaler1) { printf("Null pointer to scalers for empty class\n"); return 2; }
    if (!emptyScaler2) { printf("Null pointer to scalers for empty class\n");  return 2; }
    l0bempty+=dif(emptyScaler2->GetLOCB(),emptyScaler1->GetLOCB());
  }
  UInt_t t1 = scalers->GetScalersRecord(0         )->GetTimeStamp()->GetSeconds();
  UInt_t t2 = scalers->GetScalersRecord(nEntries-1)->GetTimeStamp()->GetSeconds();
  run_duration = dif(t2,t1);
  
  for (Int_t r=0;r<nEntries-1;r++){
    // Get SOR and EOR scaler records
    AliTriggerScalersRecord* record1 = scalers->GetScalersRecord(r);
    const AliTimeStamp*      stamp1  = record1->GetTimeStamp();
    Int_t period = stamp1->GetPeriod();
    Int_t orbit = stamp1->GetOrbit();
    printf("%5i %5i %5i %ll\n",r,period,orbit);
  }

  
  par[1] = run_duration;
  par[2] = l0b;
  if (TMath::Abs(run_duration)<1) return 3;

  // Get number of colliding bunches per orbit
  Double_t nBCsPerOrbit = -1;
  Double_t orbitRate = 11245.;
  if (refClassName.Contains("1B-ABCE-")){
    nBCsPerOrbit = lhc->GetNInteractingBunchesMeasured();
    Printf("Number of BCs from LHC data=%i",nBCsPerOrbit);
    if (nBCsPerOrbit<0 && l0bempty>0) nBCsPerOrbit = Double_t(l0bempty)/orbitRate/run_duration;
  } else {
    // Extract number of bunches per orbit
    AliTriggerClass* cl = cfg->GetTriggerClass(classId);
    AliTriggerBCMask* mask = cl->GetBCMask();
    nBCsPerOrbit = mask->GetNUnmaskedBCs();
  }
  par[3] = nBCsPerOrbit;
  
  Double_t totalBCs = orbitRate*run_duration*nBCsPerOrbit;
  if (totalBCs<1 || l0b<1) return 4;
  par[4] = -TMath::Log(1-Double_t(l0b)/totalBCs); // mu
  return 0;
}
