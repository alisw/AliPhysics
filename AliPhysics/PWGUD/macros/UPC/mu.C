// Macro to calculate pileup factor
// evgeny.kryshen@cern.ch
//
//Two cases:
//1) Reference trigger based on coincidence with beam pick-up counters to trigger beam-beam colliding bunches. 
//   Run 120072 is an example. In this case one can get total number of colliding bunches from L0b scalers for 
//   "empty" trigger corresponding to coincidence of A and C beam pick-up counters (CBEAMB-ABCE-NOPF-ALL).
//
//2) Reference trigger based on BC masks to constrain beam-beam colliding bunches (-B- in the class name). 
//   In this case total number of bunches can be computed as (run duration)*(orbit rate)*(number of B-like BCs per orbit).
//   One can switch to this option setting emptyClassName=NULL. 

UInt_t dif(UInt_t stop, UInt_t start);


Double_t mu(Int_t run=120072, char* className = "CINT1B-ABCE-NOPF-ALL", char* emptyClassName = "CBEAMB-ABCE-NOPF-ALL"){
//void mu(Int_t run=196310, char* className = "CINT7-B-NOPF-ALLNOTRD", char* emptyClassName = NULL){
  if (!TGrid::Connect("alien://")) return -1.;

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliTriggerConfiguration* cfg = (AliTriggerConfiguration*) man->Get("GRP/CTP/Config")->GetObject();

  // Get scalers
  AliTriggerRunScalers* scalers = (AliTriggerRunScalers*) man->Get("GRP/CTP/Scalers")->GetObject();
  Int_t nEntries = scalers->GetScalersRecords()->GetEntries();
  
  // Get SOR and EOR scaler records
  AliTriggerScalersRecord* record1 = scalers->GetScalersRecord(0);
  AliTriggerScalersRecord* record2 = scalers->GetScalersRecord(nEntries-1);
  
  // Extract SOR and EOR trigger counts
  Int_t classid      = cfg->GetClassIndexFromName(className);
  AliTriggerScalers* scaler1 = record1->GetTriggerScalersForClass(classid);
  AliTriggerScalers* scaler2 = record2->GetTriggerScalersForClass(classid);
  UInt_t l0b = dif(scaler2->GetLOCB(),scaler1->GetLOCB());
  
  Double_t totalBCs;
  if (emptyClassName) {
    Int_t emptyclassid = cfg->GetClassIndexFromName(emptyClassName);
    scaler1 = record1->GetTriggerScalersForClass(emptyclassid);
    scaler2 = record2->GetTriggerScalersForClass(emptyclassid);
    UInt_t l0bempty = dif(scaler2->GetLOCB(),scaler1->GetLOCB());
    totalBCs = l0bempty;
  }
  else {
    Double_t orbitRate = 11245.; // Hz
    // Extract SOR and EOR times
    AliTimeStamp* stemp1 = record1->GetTimeStamp();
    AliTimeStamp* stemp2 = record2->GetTimeStamp();
    UInt_t duration = stemp2->GetSeconds()-stemp1->GetSeconds();
    // Extract number of bunches per orbit
    AliTriggerClass* cl = cfg->GetTriggerClass(classid);
    AliTriggerBCMask* mask = cl->GetBCMask();
    Int_t nBCsPerOrbit = mask->GetNUnmaskedBCs();
    
    totalBCs = duration*orbitRate*nBCsPerOrbit;
  }
  
  Double_t mu = -TMath::Log(1-l0b/totalBCs);
  printf("l0b=%i totalBCs=%.0f mu=%f\n",l0b,totalBCs,mu);
  return mu;
}

UInt_t dif(UInt_t stop, UInt_t start){
  UInt_t d;
  if(stop >= start) d=stop-start;
  else d = stop+(0xffffffff-start)+1;
  return d;
};
