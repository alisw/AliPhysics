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

TObjArray GetClasses(Int_t run, TString ocdbStorage, TString &partition, TString &activeDetectors, Double_t& run_duration, 
    ULong64_t* LMB, ULong64_t* LMA, ULong64_t* L0B, ULong64_t* L0A, ULong64_t* L1B, ULong64_t* L1A, ULong64_t* L2B, ULong64_t* L2A,ULong64_t* class_duration){
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbStorage.Data());
  man->SetRun(run);
  // Get scalers 
  AliTriggerConfiguration* cfg = (AliTriggerConfiguration*) man->Get("GRP/CTP/Config")->GetObject();
  if (!cfg) { printf("No GRP/CTP/Config object for run %i\n",run); return TObjArray(); }
  
  partition = cfg->GetName();
  activeDetectors = cfg->GetActiveDetectors();
  
  TObjArray classes = cfg->GetClasses();
  
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*) man->Get("GRP/CTP/Scalers")->GetObject();
  if (!scalers) { printf("No GRP/CTP/Scalers object for run %i\n",run); return TObjArray(); }
  Int_t nEntries = scalers->GetScalersRecords()->GetEntriesFast();
  
  ULong64_t previousL2A[100];
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
      
      LMB[i] += dif(scaler2->GetLMCB(),scaler1->GetLMCB());
      LMA[i] += dif(scaler2->GetLMCA(),scaler1->GetLMCA());
      L0B[i] += dif(scaler2->GetLOCB(),scaler1->GetLOCB());
      L0A[i] += dif(scaler2->GetLOCA(),scaler1->GetLOCA());
      L1B[i] += dif(scaler2->GetL1CB(),scaler1->GetL1CB());
      L1A[i] += dif(scaler2->GetL1CA(),scaler1->GetL1CA());
      L2B[i] += dif(scaler2->GetL2CB(),scaler1->GetL2CB());
      L2A[i] += dif(scaler2->GetL2CA(),scaler1->GetL2CA());
      if (L2A[i]==previousL2A[i]) continue;
      previousL2A[i]=L2A[i];
      class_duration[i]+=dif(record2->GetTimeStamp()->GetSeconds(),record1->GetTimeStamp()->GetSeconds());
    }
  }
  UInt_t t1 = scalers->GetScalersRecord(0         )->GetTimeStamp()->GetSeconds();
  UInt_t t2 = scalers->GetScalersRecord(nEntries-1)->GetTimeStamp()->GetSeconds();
  run_duration = dif(t2,t1);
  return classes;
}


Int_t triggerInfo(Int_t run, TString ocdbStorage, TString &lhcPeriod, TString &lhcState, 
    Int_t &fill, Int_t &nBCsPerOrbit, Int_t &timeStart, Int_t &timeEnd){
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbStorage.Data());
  man->SetRun(run);

  AliGRPObject* grp = (AliGRPObject*) man->Get("GRP/GRP/Data")->GetObject();
  if (!grp) { printf("No GRP/GRP/Data object for run %i\n",run); return 1; }
  lhcState  = grp->GetLHCState();
  lhcPeriod = grp->GetLHCPeriod();
  timeStart = grp->GetTimeStart();
  timeEnd   = grp->GetTimeEnd();
  
  if (run==189694) return 1; // No GRP/GRP/LHCData for this run
  
  AliLHCData* lhc = (AliLHCData*) man->Get("GRP/GRP/LHCData")->GetObject();
  if (!lhc) { printf("No GRP/GRP/LHCData object for run %i\n",run); return 1; }
  fill = lhc->GetFillNumber();
  nBCsPerOrbit = lhc->GetNInteractingBunchesMeasured();

  return 0;
}
