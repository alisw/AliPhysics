void BrowseAndFillPhysicsSelectionOADB(Bool_t fill = kFALSE) {
  TString oadbfilename("../COMMON/PHYSICSSELECTION/data/physicsSelection.root");
  if (!fill) {
    TFile * f = new TFile(oadbfilename);
    new TBrowser();
    return;
  }

  AliOADBContainer * oadbContPS = new AliOADBContainer("physSel");
  AliOADBContainer * oadbContFillingScheme = new AliOADBContainer("fillScheme");
  AliOADBContainer * oadbContTriggerAnalysis = new AliOADBContainer("trigAnalysis");

  UInt_t triggerCount = 0;

  // DefaultPP (since 2012)
  AliOADBPhysicsSelection * oadbDefaultPP = new AliOADBPhysicsSelection("oadbDefaultPP");
  triggerCount = 0;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-[I|B|S]-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT5-[I|B]-NOPF-[ALL|CENT]NOTRD,C0SMB-[I|B]-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-[I|B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-[I|B]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");
  
  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT8,"+C[INT8|0TVX]-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");
  
  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kINT7inMUON,"+CINT7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");
  
  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-[B|S|SC]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt8,"+CMSH8-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt8,"+CMSL8-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt8,"+CMLL8-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt8,"+CMUL8-[B|S]-NOPF-[MUON|MUFAST|ALLNOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-[I|B|S]-NOPF-[ALL|CENT][NOTRD|],C[E|D]MC7-B-NOPF-CALOFAST","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMC8,"+C[E|D]MC8-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7[|PHL|PHM|PHH]-[I|B|S]-NOPF-[ALL|CENT][NOTRD|],CPHI7-B-NOPF-CALOFAST","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kPHI8,"+CPHI8-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-[I|B|S]-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"1");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kHighMultV0,"+CVHMV0M-B-[NOPF|SPD1|SPD2]-CENT[|NOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"VHM && V0M");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");
  
  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kHighMultSPD,"+CVHMSH2-B-[NOPF|SPD1]-CENT[|NOTRD]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"VHM && SH2");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");
  
  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kHighMultSPD,"+CSHM7-[B|S]-NOPF-ALLNOTRD","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kHighMultSPD,"+CSHM8-S-NOPF-ALLNOTRD","B",    triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+C[EMC7E|DMC7D]J[A|1|2]-[I|B|S]-NOPF-[ALL|CENT][NOTRD|],C[EMC7E|DMC7D]J[1|2]-B-NOPF-CALOFAST","B", triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+C[EMC7E|DMC7D]G[A|1|2]-[I|B|S]-NOPF-[ALL|CENT][NOTRD|],C[EMC7E|DMC7D]G[1|2]-B-NOPF-CALOFAST","B", triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+C[EMC8E|DMC8D]J[E|1|2]-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+C[EMC8E|DMC8D]G[A|1|2]-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B", triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kSPI,"+CSPI7-[B|S]-NOPF-ALLNOTRD","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kSPI,"+CSPI8-[B|S]-NOPF-[ALL|CENT][NOTRD|]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt0,"+C0MUL-[B|SA|SC]-NOPF-[MUON|MUFAST]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"1");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHJT-[B|I|S]-NOPF-[CENT|FAST]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHJT && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHSE-[B|I|S]-NOPF-[CENT|FAST]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHSE && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHQU-[B|I|S]-NOPF-[CENT|FAST]","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHQU && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CEMC7WUHQU-[B|I|S]-NOPF-CENT","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHQU && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT8WUHJT-[B|I|S]-NOPF-CENT","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHJT && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT8WUHSE-[B|I|S]-NOPF-CENT","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHSE && !IncompleteEvent");

  triggerCount++;
  oadbDefaultPP->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT8WUHQU-[B|I|S]-NOPF-CENT","B",triggerCount);
  oadbDefaultPP->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPP->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !TPCHVdip && TRDHQU && !IncompleteEvent");

  oadbContPS->AddDefaultObject(oadbDefaultPP);
  
  // DefaultPbPb (since 2015)
  AliOADBPhysicsSelection * oadbDefaultPbPb = new AliOADBPhysicsSelection("oadbDefaultPbPb");
  triggerCount = 0;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-CENT","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG && ZDCTime && !TPCHVdip");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kINT7,"+[CINT7|CV0L7]-B-NOPF-CENT","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kINT8,"+C0TVX-B-NOPF-CENT","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"T0");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"!T0BG && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"1");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kINT7inMUON,"+CINT7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUFAST","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime");
  
  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+CINT7[E|D]J[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+CINT7[E|D]G[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime && !TPCHVdip");

  triggerCount++;
  oadbDefaultPbPb->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-B-NOPF-CENTNOPMD,CINT7PH[L|M|H]-B-NOPF-CENTNOPMD,CPER7PHM-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbDefaultPbPb->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbDefaultPbPb->SetOfflineTrigger       (triggerCount,"V0A && V0C && ZDCTime && !TPCHVdip");
  
  oadbContPS->AddDefaultObject(oadbDefaultPbPb);
  
  // PbPb without ZDC
  AliOADBPhysicsSelection * oadbPbPbWithoutZDC = new AliOADBPhysicsSelection("oadbPbPbWithoutZDC");
  triggerCount = 0;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-CENT","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG");
  
  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kINT7,"+[CINT7|CV0L7]-B-NOPF-CENT","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kINT8,"+C0TVX-B-NOPF-CENT","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"T0");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"!T0BG");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kINT7inMUON,"+CINT7-B-NOPF-MUFAST","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");
  
  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUFAST","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUFAST","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");
  
  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUFAST","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");
  
  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUFAST","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");
  
  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+CINT7[E|D]J[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+CINT7[E|D]G[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbPbPbWithoutZDC->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-B-NOPF-CENTNOPMD,CINT7PH[L|M|H]-B-NOPF-CENTNOPMD,CPER7PHM-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbPbPbWithoutZDC->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbPbPbWithoutZDC->SetOfflineTrigger       (triggerCount,"V0A && V0C");
  
  oadbContPS->AppendObject(oadbPbPbWithoutZDC->Clone("oadbPbPbWithoutZDC_1"),246543,246671);
  oadbContPS->AppendObject(oadbPbPbWithoutZDC->Clone("oadbPbPbWithoutZDC_2"),244824,244889);
  oadbContPS->AppendObject(oadbPbPbWithoutZDC->Clone("oadbPbPbWithoutZDC_3"),245061,245061);
  oadbContPS->AppendObject(oadbPbPbWithoutZDC->Clone("oadbPbPbWithoutZDC_4"),245148,245148);

  
  
  // p-Pb 2016
  AliOADBPhysicsSelection * oadb_pPb_2016 = new AliOADBPhysicsSelection("oadb_pPb_2016");
  triggerCount = 0;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-CENTNOTRD","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-B-NOPF-[CENT|FAST]","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-CENTNOTRD","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");
  
  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kINT8,"+C0TVX-B-NOPF-CENT","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"T0");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");
  
  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kINT7inMUON,"+CINT7-B-NOPF-MUFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");
  
  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-B-NOPF-CENT[|NOPMD|NOTRD]","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7[|PHL|PHM|PHH]-B-NOPF-CENT[|NOPMD|NOTRD]","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-B-NOPF-UFAST","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"1");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kHighMultV0,"+CV0M7-B-NOPF-CENT","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C && V0M");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+C[EMC7E|DMC7D]J[1|2]-B-NOPF-CENT[|NOPMD|NOTRD]","B", triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+C[EMC7E|DMC7D]G[1|2]-B-NOPF-CENT[|NOPMD|NOTRD]","B", triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  triggerCount++;
  oadb_pPb_2016->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7[HJT|HSE|HQU|HNU]-T-NOPF-CENTNOPMD","B",triggerCount);
  oadb_pPb_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_pPb_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNABG");

  oadbContPS->AppendObject(oadb_pPb_2016,265304,266318);
  oadbContPS->AppendObject(oadb_pPb_2016->Clone("oadb_pPb_2016_16t"),267132,267166);


  
  // Pb-p 2016
  AliOADBPhysicsSelection * oadb_Pbp_2016 = new AliOADBPhysicsSelection("oadb_Pbp_2016");
  triggerCount = 0;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-CENTNOTRD","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-B-NOPF-[CENT|FAST]","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-CENTNOTRD","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");
  
  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kINT8,"+C0TVX-B-NOPF-CENT","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"T0");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"!T0BG && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");
  
  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kINT7inMUON,"+CINT7-B-NOPF-MUFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");
  
  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-B-NOPF-CENT[|NOPMD|NOTRD]","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7[|PHL|PHM|PHH]-B-NOPF-CENT[|NOPMD|NOTRD]","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-B-NOPF-UFAST","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"1");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG  && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kHighMultV0,"+CV0M7-B-NOPF-CENT","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C && V0M");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+C[EMC7E|DMC7D]J[1|2]-B-NOPF-CENT[|NOPMD|NOTRD]","B", triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+C[EMC7E|DMC7D]G[1|2]-B-NOPF-CENT[|NOPMD|NOTRD]","B", triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  triggerCount++;
  oadb_Pbp_2016->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7[HJT|HSE|HQU|HNU]-T-NOPF-CENTNOPMD","B",triggerCount);
  oadb_Pbp_2016->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadb_Pbp_2016->SetOfflineTrigger       (triggerCount,"V0A && V0C && !SPDClsVsTrkBG && !V0Casym && !V0C012vsTklBG && !V0MOnVsOfPileup && !SPDOnVsOfPileup && !V0PFPileup && !SPDVtxPileup && !ZNCBG");

  oadbContPS->AppendObject(oadb_Pbp_2016,266405,267131);
  
  
  // pp2010-11
  AliOADBPhysicsSelection * oadbLHCpp2010 = new AliOADBPhysicsSelection("oadbLHCpp2010");
  triggerCount=0;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1B-ABCE-NOPF-ALL,CINT1-B-NOPF-ALL[NOTRD|]","B",triggerCount);
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kINT1 | AliVEvent::kFastOnly,"+CINT1-B-NOPF-FASTNOTRD -CINT1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-[I|B]-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON,CMUS1-B-NOPF-[ALL|ALLNOTRD|MUON]","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kCMUS5,"+CMUS5-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  
  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMUSH7-B-NOPF-MUON","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMUL7-B-NOPF-MUON","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUU7-B-NOPF-[MUON|ALLNOTRD]","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kEMC7,"+CEMC7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMUS7-B-NOPF-MUON","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-[B|I]-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp");
  
  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kHighMult,"+CSH1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kEMC1,"+CEMC1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && EMCAL && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHCpp2010->AddCollisionTriggerClass(AliVEvent::kPHI1,"+CPHI1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHCpp2010->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHCpp2010->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  oadbContPS->AppendObject(oadbLHCpp2010                          ,104065,118555);
  oadbContPS->AppendObject(oadbLHCpp2010->Clone("oadbLHCpp2010_1"),118562,136377);
  oadbContPS->AppendObject(oadbLHCpp2010->Clone("oadbLHCpp2010_2"),144871,165747);
  
  // LHC10c, fill 1069 (problems with the V0 online trigger in ESD)
  AliOADBPhysicsSelection * oadbLHC10cV0Bug = new AliOADBPhysicsSelection("oadbLHC10cV0Bug");
  triggerCount=0;
  oadbLHC10cV0Bug->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1B-ABCE-NOPF-ALL","B",triggerCount);
  oadbLHC10cV0Bug->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || CTPV0A || CTPV0C");
  oadbLHC10cV0Bug->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");

  triggerCount++;
  oadbLHC10cV0Bug->AddCollisionTriggerClass(AliVEvent::kMUON,"+CMUS1B-ABCE-NOPF-MUON","B",triggerCount);
  oadbLHC10cV0Bug->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || CTPV0A || CTPV0C");
  oadbLHC10cV0Bug->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG  && !TPCLaserWarmUp");
  oadbContPS->AppendObject(oadbLHC10cV0Bug, 118556,118561);

  // LHC10h
  AliOADBPhysicsSelection * oadbLHC10h = new AliOADBPhysicsSelection("oadbLHC10h");
  triggerCount=0;
  oadbLHC10h->AddCollisionTriggerClass(AliVEvent::kMB,"+CMBAC-B-NOPF-ALL","B",triggerCount);
  oadbLHC10h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC10h->SetOfflineTrigger       (triggerCount,"V0A && V0C && SPDGFOL1 > 1 && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC10h->AddCollisionTriggerClass(AliVEvent::kMB,"+CMBS2A-B-NOPF-ALL","B",triggerCount);
  oadbLHC10h->SetHardwareTrigger      (triggerCount,"SPDGFOL1 > 1 && V0A");
  oadbLHC10h->SetOfflineTrigger       (triggerCount,"V0A && V0C && SPDGFOL1 > 1 && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC10h->AddCollisionTriggerClass(AliVEvent::kMB,"+CMBS2C-B-NOPF-ALL","B",triggerCount);
  oadbLHC10h->SetHardwareTrigger      (triggerCount,"SPDGFOL1 > 1 && V0C");
  oadbLHC10h->SetOfflineTrigger       (triggerCount,"V0A && V0C && SPDGFOL1 > 1 && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC10h->AddCollisionTriggerClass(AliVEvent::kMB,"+CMBACS2-B-NOPF-ALL[NOTRD|]","B",triggerCount);
  oadbLHC10h->SetHardwareTrigger      (triggerCount,"V0A && V0C && SPDGFOL1 > 1");
  oadbLHC10h->SetOfflineTrigger       (triggerCount,"V0A && V0C && SPDGFOL1 > 1 && !TPCLaserWarmUp && ZDCTime");
  
  triggerCount++;
  oadbLHC10h->AddCollisionTriggerClass(AliVEvent::kHighMult,"+C0SMH-B-NOPF-ALL[NOTRD|]","B",triggerCount);
  oadbLHC10h->SetHardwareTrigger      (triggerCount,"SPDGFO >= 100");
  oadbLHC10h->SetOfflineTrigger       (triggerCount,"SPDGFO >= 100 && !V0ABG && !V0CBG && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC10h,136851,139517);

  
  // LHC11h
  AliOADBPhysicsSelection * oadbLHC11h = new AliOADBPhysicsSelection("oadbDefaultlhc11h");
  triggerCount=0;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kMB,"+CPBI[1|2_B1]-B-[NOPF|PF]-ALLNOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kCentral,"+C[SEMI|CENT|VLN|VHN][|_B2|_R1|_R2]-B-[NOPF|PF]-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C && Central");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");
  
  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kSemiCentral,"+C[SEMI|CENT|VLN|VHN][|_B2|_R1|_R2]-B-[NOPF|PF]-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C && SemiCentral && !Central");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+CPBI2EJE-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+CPBI2EGA-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kMUSPB,"+CPBI1MSL-B-NOPF-MUON","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kMUSHPB,"+CPBI1MSH-B-NOPF-MUON","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kMuonUnlikePB,"+CPBI1MUL-B-NOPF-MUON","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kMuonLikePB,"+CPBI1MLL-B-NOPF-MUON","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  triggerCount++;
  oadbLHC11h->AddCollisionTriggerClass(AliVEvent::kPHOSPb,"+CPBI2PHS-B-NOPF-CENTNOTRD","B",triggerCount);
  oadbLHC11h->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC11h->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCLaserWarmUp && ZDCTime");

  oadbContPS->AppendObject(oadbLHC11h, 166529, 170593);
  
  // LHC15m: quiet Pb-Pb without ZDC selection
  AliOADBPhysicsSelection * oadbLHC15m = new AliOADBPhysicsSelection("oadbLHC15m");
  triggerCount=0;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-CENTNOTRD","B", triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !TPCHVdip");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-CENT","B", triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG && !TPCHVdip");
  
  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-B-NOPF-CENT","B", triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCHVdip");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kINT8,"+C0TVX-B-NOPF-CENT","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"T0");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"!T0BG && !TPCHVdip");
  
  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUFAST","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUFAST","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUFAST","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUFAST","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kEMC7,"+C[E|D]MC7-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCHVdip");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCHVdip");
  
  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+C[EMC7E|DMC7D]J[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCHVdip");

  triggerCount++;
  oadbLHC15m->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+C[EMC7E|DMC7D]G[1|2]-B-NOPF-CENTNOPMD","B",triggerCount);
  oadbLHC15m->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC15m->SetOfflineTrigger       (triggerCount,"V0A && V0C && !TPCHVdip");
  
  oadbContPS->AppendObject(oadbLHC15m, 243890, 243984);

  
  // p-Pb (Pb going in A direction)
  AliOADBPhysicsSelection * oadbLHC13b = new AliOADBPhysicsSelection("oadbLHC13b");
  triggerCount=0;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kINT5,"+CINT5-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A || V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"(V0A || V0C) && !V0ABG && !V0CBG && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-[I|B]-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kSPI7,"+CSPI7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kEMC7,"+CEMC7-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+CEMC7EG[1|2]-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+CEMC7EJ[1|2]-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-[MUON|ALLNOTRD]","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"1");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG && !ZNABG && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHJT-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHJT-B-NOPF-FAST -CINT7WUHJT-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && TRDHJT && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHSE-B-NOPF-[ALL|CENT]", "B",triggerCount);
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHSE-B-NOPF-FAST -CINT7WUHSE-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && TRDHSE && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHQU-B-NOPF-[ALL|CENT]", "B",   triggerCount);
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHQU-B-NOPF-FAST -CINT7WUHQU-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && TRDHQU && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kTRD,"+CEMC7WUHEE-B-NOPF-CENT", "B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNABG && TRDHEE && !TPCHVdip && !IncompleteEvent");

  triggerCount++;
  oadbLHC13b->AddCollisionTriggerClass(AliVEvent::kHighMult,"+CSHM7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13b->SetHardwareTrigger      (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C");
  oadbLHC13b->SetOfflineTrigger       (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C && !ZNABG && !TPCHVdip && !IncompleteEvent");

  oadbContPS->AppendObject(oadbLHC13b                       ,194713,196345);
  oadbContPS->AppendObject(oadbLHC13b->Clone("oadbLHC13b_2"),188124,188374);

  // Pb-p (Pb going in C direction)
  AliOADBPhysicsSelection * oadbLHC13f = new AliOADBPhysicsSelection("oadbLHC13f");
  triggerCount=0;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kINT1,"+CINT1-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"SPDGFO >= 1 || V0A || V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"(SPDGFO >= 1 || V0A || V0C) && !V0ABG && !V0CBG && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kINT7,"+CINT7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kEMC7,"+CEMC7-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kEMCEGA,"+CEMC7EG[1|2]-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kEMCEJE,"+CEMC7EJ[1|2]-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kPHI7,"+CPHI7-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kMuonSingleLowPt7,"+CMSL7-B-NOPF-MUON,CMSL7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kMuonSingleHighPt7,"+CMSH7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kMuonLikeLowPt7,"+CMLL7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kMuonUnlikeLowPt7,"+CMUL7-B-NOPF-MUON","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kZED,"+C1ZED-B-NOPF-[ALL|CENT]NOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"1");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"(ZDCTDCA || ZDCTDCC) && !V0ABG && !V0CBG && !ZNCBG && !SPDClsVsTrkBG && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHJT-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHJT-B-NOPF-FAST -CINT7WUHJT-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && TRDHJT && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHSE-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHSE-B-NOPF-FAST -CINT7WUHSE-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && TRDHSE && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD,"+CINT7WUHQU-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD | AliVEvent::kFastOnly,"+CINT7WUHQU-B-NOPF-FAST -CINT7WUHQU-B-NOPF-[ALL|CENT]","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && TRDHQU && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kTRD,"+CEMC7WUHEE-B-NOPF-CENT","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"V0A && V0C && !ZNCBG && TRDHEE && !TPCHVdip");

  triggerCount++;
  oadbLHC13f->AddCollisionTriggerClass(AliVEvent::kHighMult,"+CSHM7-B-NOPF-ALLNOTRD","B",triggerCount);
  oadbLHC13f->SetHardwareTrigger      (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C");
  oadbLHC13f->SetOfflineTrigger       (triggerCount,"SPDGFOL1 >= 10 && V0A && V0C && !ZNCBG && !TPCHVdip");

  oadbContPS->AppendObject(oadbLHC13f, 196346, 197411);

  // filling schemes
  // Defaults
  AliOADBFillingScheme * fsDefault = new AliOADBFillingScheme("Default");
  fsDefault->SetFillingSchemeName("Default");
  fsDefault->SetBXIDs("B","");
  fsDefault->SetBXIDs("A","");
  fsDefault->SetBXIDs("AC","");
  fsDefault->SetBXIDs("ACE","");
  fsDefault->SetBXIDs("C","");
  fsDefault->SetBXIDs("E","");
  oadbContFillingScheme->AddDefaultObject(fsDefault);

  AliOADBFillingScheme * fs4by4a = new AliOADBFillingScheme("4x4a");
  fs4by4a->SetFillingSchemeName("4x4a");
  fs4by4a->SetBXIDs("B"," #2128 #3019");
  fs4by4a->SetBXIDs("A"," #346 #3465" );
  fs4by4a->SetBXIDs("AC","");
  fs4by4a->SetBXIDs("C"," #1234 #1680");
  fs4by4a->SetBXIDs("E"," #790");
  oadbContFillingScheme->AppendObject(fs4by4a, 104065, 104160);

  AliOADBFillingScheme * fs4by4astar = new AliOADBFillingScheme("4x4a*");
  fs4by4astar->SetFillingSchemeName("4x4a*");
  fs4by4astar->SetBXIDs("B"," #2000 #2891");
  fs4by4astar->SetBXIDs("A"," #218 #3337" );
  fs4by4astar->SetBXIDs("AC","");
  fs4by4astar->SetBXIDs("C"," #1106 #1552");
  fs4by4astar->SetBXIDs("E"," #790");
  oadbContFillingScheme->AppendObject(fs4by4astar, 104315, 104321);

  AliOADBFillingScheme * fs4by4b = new AliOADBFillingScheme("4x4b");
  fs4by4b->SetFillingSchemeName("4x4b");
  fs4by4b->SetBXIDs("B"," #2228 #3119");
  fs4by4b->SetBXIDs("A"," #2554 #446" );
  fs4by4b->SetBXIDs("AC","");
  fs4by4b->SetBXIDs("C"," #1334 #769");
  fs4by4b->SetBXIDs("E"," #790");
  oadbContFillingScheme->AppendObject(fs4by4b, 104792, 104803);

  AliOADBFillingScheme * fs4by4c = new AliOADBFillingScheme("4x4c");
  fs4by4c->SetFillingSchemeName("4x4c");
  fs4by4c->SetBXIDs("B"," #3119 #769");
  fs4by4c->SetBXIDs("A"," #2554 #446" );
  fs4by4c->SetBXIDs("AC"," ");
  fs4by4c->SetBXIDs("C"," #1334 #2228");
  fs4by4c->SetBXIDs("E"," #790");
  oadbContFillingScheme->AppendObject(fs4by4c, 104824, 104892);

  AliOADBFillingScheme * fs16by16a = new AliOADBFillingScheme("16x16a");
  fs16by16a->SetFillingSchemeName("16x16a");
  fs16by16a->SetBXIDs("B"," #1337 #1418 #2228 #2309 #3119 #3200 #446 #527");
  fs16by16a->SetBXIDs("A"," #1580  #1742  #1904  #2066  #2630  #2792  #2954  #3362" );
  fs16by16a->SetBXIDs("AC","");
  fs16by16a->SetBXIDs("C"," #845  #1007  #1169   #1577 #3359 #3521 #119  #281");
  fs16by16a->SetBXIDs("E"," #790");
  oadbContFillingScheme->AppendObject(fs16by16a, 105143, 105160);
  oadbContFillingScheme->AppendObject(fs4by4c->Clone(), 105256, 105268);

  AliOADBFillingScheme * fsSingle_2b_1_1_1 = new AliOADBFillingScheme("Single_2b_1_1_1");
  fsSingle_2b_1_1_1->SetFillingSchemeName("Single_2b_1_1_1");
  fsSingle_2b_1_1_1->SetBXIDs("B"," #346");
  fsSingle_2b_1_1_1->SetBXIDs("A"," #2131" );
  fsSingle_2b_1_1_1->SetBXIDs("AC","");
  fsSingle_2b_1_1_1->SetBXIDs("C"," #3019");
  fsSingle_2b_1_1_1->SetBXIDs("E"," #1238");
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1, 114786, 116684);

  AliOADBFillingScheme * fsSingle_3b_2_2_2 = new AliOADBFillingScheme("Single_3b_2_2_2");
  fsSingle_3b_2_2_2->SetFillingSchemeName("Single_3b_2_2_2");
  fsSingle_3b_2_2_2->SetBXIDs("B","   #346  #1240");
  fsSingle_3b_2_2_2->SetBXIDs("A","  #2131" );
  fsSingle_3b_2_2_2->SetBXIDs("AC","");
  fsSingle_3b_2_2_2->SetBXIDs("C"," #3019 ");
  fsSingle_3b_2_2_2->SetBXIDs("E"," #1238");
  oadbContFillingScheme->AppendObject(fsSingle_3b_2_2_2, 117048, 117120);
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1->Clone(), 117220, 118555);

  AliOADBFillingScheme * fsSingle_2b_1_1_1_mis = new AliOADBFillingScheme("Single_2b_1_1_1 - 1 misaligned bx");
  fsSingle_2b_1_1_1_mis->SetFillingSchemeName("Single_2b_1_1_1 - 1 misaligned BX");
  fsSingle_2b_1_1_1_mis->SetBXIDs("B"," #345");
  fsSingle_2b_1_1_1_mis->SetBXIDs("A"," #2130" );
  fsSingle_2b_1_1_1_mis->SetBXIDs("AC","");
  fsSingle_2b_1_1_1_mis->SetBXIDs("C"," #3018");
  fsSingle_2b_1_1_1_mis->SetBXIDs("E"," #1238");
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1_mis, 118556, 118783);
  oadbContFillingScheme->AppendObject(fsSingle_2b_1_1_1->Clone(), 118784, 119163);

  AliOADBFillingScheme * fsSingle_4b_2_2_2 = new AliOADBFillingScheme("Single_4b_2_2_2");
  fsSingle_4b_2_2_2->SetFillingSchemeName("Single_4b_2_2_2");
  fsSingle_4b_2_2_2->SetBXIDs("B","   #669  #3019 ");
  fsSingle_4b_2_2_2->SetBXIDs("A","  #346  #2454 " );
  fsSingle_4b_2_2_2->SetBXIDs("AC","");
  fsSingle_4b_2_2_2->SetBXIDs("C","  #1234  #2128 ");
  fsSingle_4b_2_2_2->SetBXIDs("E"," #1681 #3463");
  oadbContFillingScheme->AppendObject(fsSingle_4b_2_2_2, 119837, 119862);

  AliOADBFillingScheme * fsSingle_6b_3_3_3 = new AliOADBFillingScheme("Single_6b_3_3_3");
  fsSingle_6b_3_3_3->SetFillingSchemeName("Single_6b_3_3_3");
  fsSingle_6b_3_3_3->SetBXIDs("B","   #346  #546  #746 ");
  fsSingle_6b_3_3_3->SetBXIDs("A","  #2131  #2331  #2531 " );
  fsSingle_6b_3_3_3->SetBXIDs("AC","");
  fsSingle_6b_3_3_3->SetBXIDs("C"," #3019  #3219  #3419");
  fsSingle_6b_3_3_3->SetBXIDs("E"," #1296 #1670");
  oadbContFillingScheme->AppendObject(fsSingle_6b_3_3_3, 119902, 120691);

  AliOADBFillingScheme * fsSingle_13b_8_8_8 = new AliOADBFillingScheme("Single_13b_8_8_8");
  fsSingle_13b_8_8_8->SetFillingSchemeName("Single_13b_8_8_8");
  fsSingle_13b_8_8_8->SetBXIDs("B","  #346  #446  #546  #646  #1240  #1340  #1440  #1540");
  fsSingle_13b_8_8_8->SetBXIDs("A","  #946  #2131  #2231  #2331  #2431 " );
  fsSingle_13b_8_8_8->SetBXIDs("AC","");
  fsSingle_13b_8_8_8->SetBXIDs("C"," #3019  #3119  #3219  #3319  #3519 ");
  fsSingle_13b_8_8_8->SetBXIDs("E"," #1835 #2726");
  oadbContFillingScheme->AppendObject(fsSingle_13b_8_8_8, 120741, 122375);

  AliOADBFillingScheme * fs125n_48b_36_16_36 = new AliOADBFillingScheme("125n_48b_36_16_36");
  fs125n_48b_36_16_36->SetFillingSchemeName("125n_48b_36_16_36");
  fs125n_48b_36_16_36->SetBXIDs("B","   #346  #396  #446  #496  #546  #596  #646  #696  #1240  #1290  #1340  #1390  #1440  #1490  #1540  #1590 ");
  fs125n_48b_36_16_36->SetBXIDs("A"," #755  #805  #855  #905  #955  #1005  #1799  #1849  #1899  #2131  #2181  #2231  #2281  #2331  #2381  #2431  #2481  #2531  #2581  #2631  #2846  #3016  #3066  #3116  #3166  #3216  #3266  #3316  #3366  #3425  #3475  #3525 " );
  fs125n_48b_36_16_36->SetBXIDs("AC","");
  fs125n_48b_36_16_36->SetBXIDs("C","  #3019  #3069  #3119  #3169  #3219  #3269  #3319  #3369  #14  #64  #114  #746  #796  #846  #908  #958  #1008  #1640  #1690  #1740  #2055  #2125  #2175  #2225  #2275  #2325  #2375  #2425  #2475  #2534  #2584  #2634 ");
  fs125n_48b_36_16_36->SetBXIDs("E","");
  oadbContFillingScheme->AppendObject(fs125n_48b_36_16_36, 130148, 130375);

  AliOADBFillingScheme * fs1000ns_50b_35_14_35 = new AliOADBFillingScheme("1000ns_50b_35_14_35");
  fs1000ns_50b_35_14_35->SetFillingSchemeName("1000ns_50b_35_14_35");
  fs1000ns_50b_35_14_35->SetBXIDs("B","  #346  #386  #426  #466  #506  #546  #586  #1240  #1280  #1320  #1360  #1400  #1440  #1480 ");
  fs1000ns_50b_35_14_35->SetBXIDs("A","  #626  #666  #706  #746  #786  #826  #866  #1520  #1560  #1600  #1640  #1680  #1720  #1760  #2076  #2131  #2171  #2211  #2251  #2291  #2331  #2371  #2414  #2454  #2494  #2534  #2574  #2614  #2654  #2694  #2734  #2774  #2814 " ); //#2854  #2894  #2934 not present in this run
  fs1000ns_50b_35_14_35->SetBXIDs("AC","");
  fs1000ns_50b_35_14_35->SetBXIDs("C"," #3019  #3059  #3099  #3139  #3179  #3219  #3259  #3299  #3339  #3379  #3419  #3459  #3499  #3539  #115  #629  #669  #709  #749  #789  #829  #869  #909  #949  #989  #1029  #1069  #1109  #1149  #1523  #1563  #1603  #1643 "); //#1683  #1723  #1763 not present in this run
  fs1000ns_50b_35_14_35->SetBXIDs("E","");
  oadbContFillingScheme->AppendObject(fs1000ns_50b_35_14_35, 130601, 130640);


  AliOADBFillingScheme * fs200ns_192Pb_216p_9inj24bpi = new AliOADBFillingScheme("200ns_192Pb_216p_9inj24bpi");
  fs200ns_192Pb_216p_9inj24bpi->SetFillingSchemeName("200ns_192Pb_216p_9inj24bpi");
  fs200ns_192Pb_216p_9inj24bpi->SetBXIDs("B","  #346  #354  #363  #371  #380  #388  #397  #405  #414  #422  #431  #439  #448  #456  #465  #473  #482  #490  #499  #507  #516  #524  #533  #541  #578  #586  #595  #603  #612  #620  #629  #637  #646  #654  #663  #671  #680  #688  #697  #705  #714  #722  #731  #739  #748  #756  #765  #773  #810  #818  #827  #835  #844  #852  #861  #869  #878  #886  #895  #903  #912  #920  #929  #937  #946  #954  #963  #971  #980  #988  #997  #1005  #1237  #1245  #1254  #1262  #1271  #1279  #1288  #1296  #1305  #1313  #1322  #1330  #1339  #1347  #1356  #1364  #1373  #1381  #1390  #1398  #1407  #1415  #1424  #1432  #1469  #1477  #1486  #1494  #1503  #1511  #1520  #1528  #1537  #1545  #1554  #1562  #1571  #1579  #1588  #1596  #1605  #1613  #1622  #1630  #1639  #1647  #1656  #1664  #1701  #1709  #1718  #1726  #1735  #1743  #1752  #1760  #1769  #1777  #1786  #1794  #1803  #1811  #1820  #1828  #1837  #1845  #1854  #1862  #1871  #1879  #1888  #1896  #2128  #2136  #2145  #2153  #2162  #2170  #2179  #2187  #2196  #2204  #2213  #2221  #2230  #2238  #2247  #2255  #2264  #2272  #2281  #2289  #2298  #2306  #2315  #2323  #3019  #3027  #3036  #3044  #3053  #3061  #3070  #3078  #3087  #3095  #3104  #3112  #3121  #3129  #3138  #3146  #3155  #3163  #3172  #3180  #3189  #3197  #3206  #3214 ");
  fs200ns_192Pb_216p_9inj24bpi->SetBXIDs("A","");
  fs200ns_192Pb_216p_9inj24bpi->SetBXIDs("C","  #1  #9  #18  #26  #35  #43  #52  #60  #69  #77  #86  #94  #103  #111  #120  #128  #3497  #3505  #3514  #3522  #3531  #3539  #3548  #3556 ");
  fs200ns_192Pb_216p_9inj24bpi->SetBXIDs("E","  #176  #261 ");
  oadbContFillingScheme->AppendObject(fs200ns_192Pb_216p_9inj24bpi, 196474, 196474);

  // Trigger analysis
  AliOADBTriggerAnalysis * oadbTrigAnalysis = new AliOADBTriggerAnalysis("Default");
  oadbTrigAnalysis->SetZDCCorrParameters(0., 0., 2., 2.);
  oadbContTriggerAnalysis->AddDefaultObject(oadbTrigAnalysis);

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC1 = new AliOADBTriggerAnalysis("ZDCCut1");
  oadbTrigAnalysisZDC1->SetZDCCorrParameters(-66.9, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC1->Clone("ZDCCut1_pass1"), 136851, 137848,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC1->Clone("ZDCCut1_pass2"), 136851, 137848,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC2 = new AliOADBTriggerAnalysis("ZDCCut2");
  oadbTrigAnalysisZDC2->SetZDCCorrParameters(-66.2, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC2->Clone("ZDCCut2_pass1"), 138125, 138275,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC2->Clone("ZDCCut2_pass2"), 138125, 138275,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC3 = new AliOADBTriggerAnalysis("ZDCCut3");
  oadbTrigAnalysisZDC3->SetZDCCorrParameters(-65.4, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC3->Clone("ZDCCut3_pass1"), 138359, 138469,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC3->Clone("ZDCCut3_pass2"), 138359, 138469,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC4 = new AliOADBTriggerAnalysis("ZDCCut4");
  oadbTrigAnalysisZDC4->SetZDCCorrParameters(-67.7, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC4->Clone("ZDCCut4_pass1"), 138533, 138742,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC4->Clone("ZDCCut4_pass2"), 138533, 138742,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC5 = new AliOADBTriggerAnalysis("ZDCCut5");
  oadbTrigAnalysisZDC5->SetZDCCorrParameters(-67.2, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC5->Clone("ZDCCut5_pass1"), 138795, 138872,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC5->Clone("ZDCCut5_pass2"), 138795, 138872,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC6 = new AliOADBTriggerAnalysis("ZDCCut6");
  oadbTrigAnalysisZDC6->SetZDCCorrParameters(-65.6, -2.1, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC6->Clone("ZDCCut6_pass1"), 138924, 139517,"pass1");
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC6->Clone("ZDCCut6_pass2"), 138924, 139517,"pass2");

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC7 = new AliOADBTriggerAnalysis("ZDCCut7");
  oadbTrigAnalysisZDC7->SetZNCorrParameters(2.0,100,5.,100);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC7, 194713, 196345);

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC8 = new AliOADBTriggerAnalysis("ZDCCut8");
  oadbTrigAnalysisZDC8->SetZNCorrParameters(5.0,100,2.0,100);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC8, 196346, 197411);

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDCCut_pPb_2016 = new AliOADBTriggerAnalysis("ZDCCut_pPb_2016");
  oadbTrigAnalysisZDCCut_pPb_2016->SetZNCorrParameters(2.0,100,5.,100);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDCCut_pPb_2016, 265304,266318);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDCCut_pPb_2016->Clone("ZDCCut_pPb_2016_16t"),267132,267166);

  AliOADBTriggerAnalysis * oadbTrigAnalysisZDCCut_Pbp_2016 = new AliOADBTriggerAnalysis("ZDCCut_Pbp_2016");
  oadbTrigAnalysisZDCCut_Pbp_2016->SetZNCorrParameters(5.0,100,2.0,100);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDCCut_Pbp_2016, 266405,267131);
  
  AliOADBTriggerAnalysis * oadbTrigAnalysisZDC9 = new AliOADBTriggerAnalysis("ZDCCut9");
  oadbTrigAnalysisZDC9->SetZDCCorrParameters(-2.1, 0, 4*0.58, 4*0.5);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisZDC9, 136851, 139517);

  AliOADBTriggerAnalysis * oadbTrigAnalysisLHC11h = new AliOADBTriggerAnalysis("lhc11h");
  oadbTrigAnalysisLHC11h->SetZDCCorrParameters(0.5, 0, 4*0.7, 4*0.7);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC11h, 166529, 170593);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15f1 = new AliOADBTriggerAnalysis("lhc15f_isolated_bunches");
  oadbTrigAnalysisLHC15f1->SetV0MOnVsOfA(0.);
  oadbTrigAnalysisLHC15f1->SetV0MOnVsOfB(0.);
  oadbTrigAnalysisLHC15f1->SetSPDOnVsOfA(0.);
  oadbTrigAnalysisLHC15f1->SetSPDOnVsOfB(0.);
  oadbTrigAnalysisLHC15f1->SetVIRBBAflags(33);
  oadbTrigAnalysisLHC15f1->SetVIRBBCflags(33);
  oadbTrigAnalysisLHC15f1->SetV0CasymA(0);
  oadbTrigAnalysisLHC15f1->SetV0CasymB(0);
  oadbTrigAnalysisLHC15f1->SetNBCsPast(0);
  oadbTrigAnalysisLHC15f1->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15f1                                   ,225000,225719);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15f1->Clone("lhc15f_isolated_bunches2"),226062,226500);

  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15f2 = new AliOADBTriggerAnalysis("lhc15f_50ns_trains1");
  oadbTrigAnalysisLHC15f2->SetV0MOnVsOfA(-372.579114);
  oadbTrigAnalysisLHC15f2->SetV0MOnVsOfB(9.415265);
  oadbTrigAnalysisLHC15f2->SetSPDOnVsOfA(-6.65857);
  oadbTrigAnalysisLHC15f2->SetSPDOnVsOfB(0.546801);
  oadbTrigAnalysisLHC15f2->SetNBCsPast(0);
  oadbTrigAnalysisLHC15f2->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15f2,225753,225768);

  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15f3 = new AliOADBTriggerAnalysis("lhc15f_50ns_trains2_missing_V0C3");
  oadbTrigAnalysisLHC15f3->SetV0MOnVsOfA(-372.579114);
  oadbTrigAnalysisLHC15f3->SetV0MOnVsOfB(9.415265);
  oadbTrigAnalysisLHC15f3->SetSPDOnVsOfA(-6.65857);
  oadbTrigAnalysisLHC15f3->SetSPDOnVsOfB(0.546801);
  oadbTrigAnalysisLHC15f3->SetV0CasymA(0);
  oadbTrigAnalysisLHC15f3->SetV0CasymB(0);
  oadbTrigAnalysisLHC15f3->SetNBCsPast(0);
  oadbTrigAnalysisLHC15f3->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15f3,226530,226606);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15h1 = new AliOADBTriggerAnalysis("lhc15h");
  oadbTrigAnalysisLHC15h1->SetV0MOnVsOfA(-245.12);
  oadbTrigAnalysisLHC15h1->SetV0MOnVsOfB(6.86754);
  oadbTrigAnalysisLHC15h1->SetSPDOnVsOfA(-6.65857);
  oadbTrigAnalysisLHC15h1->SetSPDOnVsOfB(0.546801);
  oadbTrigAnalysisLHC15h1->SetNBCsPast(0);
  oadbTrigAnalysisLHC15h1->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15h1,232914,233911);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15h2 = new AliOADBTriggerAnalysis("lhc15h_isolated_bunches");
  oadbTrigAnalysisLHC15h2->SetV0MOnVsOfA(0);
  oadbTrigAnalysisLHC15h2->SetV0MOnVsOfB(0);
  oadbTrigAnalysisLHC15h2->SetSPDOnVsOfA(0);
  oadbTrigAnalysisLHC15h2->SetSPDOnVsOfB(0);
  oadbTrigAnalysisLHC15h2->SetVIRBBAflags(33);
  oadbTrigAnalysisLHC15h2->SetVIRBBCflags(33);
  oadbTrigAnalysisLHC15h2->SetNBCsPast(0);
  oadbTrigAnalysisLHC15h2->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15h2,233912,234050);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15i = new AliOADBTriggerAnalysis("lhc15i");
  oadbTrigAnalysisLHC15i->SetV0MOnVsOfA(-223.155660);
  oadbTrigAnalysisLHC15i->SetV0MOnVsOfB(7.117266);
  oadbTrigAnalysisLHC15i->SetSPDOnVsOfA(-6.218793);
  oadbTrigAnalysisLHC15i->SetSPDOnVsOfB(0.543201);
  oadbTrigAnalysisLHC15i->SetNBCsPast(0);
  oadbTrigAnalysisLHC15i->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15i,235193,236866);

  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15j = new AliOADBTriggerAnalysis("lhc15j");
  oadbTrigAnalysisLHC15j->SetV0MOnVsOfA(-222.631866);
  oadbTrigAnalysisLHC15j->SetV0MOnVsOfB(7.431432);
  oadbTrigAnalysisLHC15j->SetSPDOnVsOfA(-6.610850);
  oadbTrigAnalysisLHC15j->SetSPDOnVsOfB(0.587165);
  oadbTrigAnalysisLHC15j->SetNBCsPast(0);
  oadbTrigAnalysisLHC15j->SetNBCsFuture(11);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15j,236892,238621);

  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15l = new AliOADBTriggerAnalysis("lhc15l");
  oadbTrigAnalysisLHC15l->SetV0MOnVsOfA(-198.639921);
  oadbTrigAnalysisLHC15l->SetV0MOnVsOfB(7.454714);
  oadbTrigAnalysisLHC15l->SetSPDOnVsOfA(-5.018572);
  oadbTrigAnalysisLHC15l->SetSPDOnVsOfB(0.585245);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15l,239188,241544);
  

  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15n = new AliOADBTriggerAnalysis("lhc15n");
  oadbTrigAnalysisLHC15n->SetV0MOnVsOfA(-336.279729);
  oadbTrigAnalysisLHC15n->SetV0MOnVsOfB(10.694535);
  oadbTrigAnalysisLHC15n->SetSPDOnVsOfA(-4.144493);
  oadbTrigAnalysisLHC15n->SetSPDOnVsOfB(0.851104); 
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15n,244340,244628);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC15o = new AliOADBTriggerAnalysis("lhc15o1");
  oadbTrigAnalysisLHC15o->SetV0MOnVsOfA(0.);
  oadbTrigAnalysisLHC15o->SetV0MOnVsOfB(0.);
  oadbTrigAnalysisLHC15o->SetSPDOnVsOfA(0.);
  oadbTrigAnalysisLHC15o->SetSPDOnVsOfB(0.);
  oadbTrigAnalysisLHC15o->SetVIRBBAflags(33);
  oadbTrigAnalysisLHC15o->SetVIRBBCflags(33);
  oadbTrigAnalysisLHC15o->SetV0CasymA(0);
  oadbTrigAnalysisLHC15o->SetV0CasymB(0);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15o,244824,245725);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15o->Clone("lhc15o2"),245794,246994);
  
  AliOADBTriggerAnalysis * oadbTrigAnalysisLHC15o1 = oadbTrigAnalysisLHC15o->Clone("lhc15o_common_zna_tdc");
  oadbTrigAnalysisLHC15o1->SetZDCCorrParameters(-123.1, 123.1, 2., 2.);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC15o1, 245726, 245793);
  
  AliOADBTriggerAnalysis* oadbTrigAnalysisLHC16do = new AliOADBTriggerAnalysis("lhc16do");
  oadbTrigAnalysisLHC16do->SetV0MOnVsOfA(-65.42);
  oadbTrigAnalysisLHC16do->SetV0MOnVsOfB(7.43);
  oadbContTriggerAnalysis->AppendObject(oadbTrigAnalysisLHC16do, 252235, 264035);
  
  oadbTrigAnalysisZDC1->Print();
  oadbTrigAnalysisZDC2->Print();
  oadbTrigAnalysisZDC3->Print();
  oadbTrigAnalysisZDC4->Print();
  oadbTrigAnalysisZDC5->Print();
  oadbTrigAnalysisZDC6->Print();
  oadbTrigAnalysisZDC7->Print();
  oadbTrigAnalysisZDC8->Print();
  oadbTrigAnalysisZDC9->Print();
  oadbTrigAnalysisLHC11h->Print();

  oadbContPS->WriteToFile(oadbfilename.Data());
  oadbContFillingScheme->WriteToFile(oadbfilename.Data());
  oadbContTriggerAnalysis->WriteToFile(oadbfilename.Data());

  TFile * fopen = new TFile (oadbfilename);
  new TBrowser;
}
