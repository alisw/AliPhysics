AliAnalysisTask *AddTaskTender(Bool_t useV0=kFALSE, 
                               Bool_t useTPC=kTRUE,
                               Bool_t useTOF=kTRUE,
                               Bool_t useTRD=kTRUE,
                               Bool_t usePID=kTRUE,
                               Bool_t useVTX=kTRUE,
                               Bool_t useT0=kTRUE,
                               Bool_t useEmc=kFALSE,
                               Bool_t usePtFix=kFALSE)
{
  if (!(useV0 | useTPC | useTOF | useTRD | usePID | useVTX | | useT0 | useEmc | usePtFix)) {
     ::Error("AddTaskTender", "No supply added to tender, so tender not created");
     return 0;
  }   
  //get the current analysis manager
  Bool_t checkEvtSelection = useV0;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTender", "No analysis manager found.");
    return 0;
  }
  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliESDInputHandler::Class())) {
    ::Error("AddTask_tender_Tender","The analysis tender only works with ESD input!");
    return 0;
  }
  
  //========= Add tender to the ANALYSIS manager and set default storage =====
  AliTender *tender=new AliTender("AnalysisTender");
  tender->SetCheckEventSelection(checkEvtSelection);
  tender->SetDefaultCDBStorage("raw://");
  mgr->AddTask(tender);
  
  //check that that tender is the first task after the pid response
  TString firstName(mgr->GetTasks()->First()->GetName());
  Bool_t isSecond=(mgr->GetTasks()->At(1) == (TObject*)tender);

  if (! (firstName=="PIDResponseTask" && isSecond )){
    Fatal("AddTaskTender","When using the tender the first task needs to be the PIDResponse and the tender the second task!!!");
    return NULL;
  }
  
  //========= Attach VZERO supply ======
  if (useV0) {
     AliVZEROTenderSupply *vzeroSupply=new AliVZEROTenderSupply("VZEROtender");
     vzeroSupply->SetDebug(kFALSE);
     tender->AddSupply(vzeroSupply);
  }   

 
  //========= Attach TPC supply ======
  if (useTPC) {
     AliTPCTenderSupply *tpcSupply=new AliTPCTenderSupply("TPCtender");
     tpcSupply->SetDebugLevel(2);
     //tpcSupply->SetMip(50.);
     tender->AddSupply(tpcSupply);
  }   

  //========= Attach track 1/pt correction supply ======
  if (usePtFix) {
     AliTrackFixTenderSupply *trfixSupply=new AliTrackFixTenderSupply("PTInvFix");
     //trfixSupply->SetDebugLevel(2);
     tender->AddSupply(trfixSupply);
  }   

  //========= Attach T0 supply ======
  if (useT0) {
    AliT0TenderSupply *t0Tender = new AliT0TenderSupply("T0tender");
    t0Tender ->SetPass4LHC11aCorrection(kTRUE);
    tender->AddSupply(t0Tender);
  }   

  //========= Attach TOF supply ======
  if (useTOF) {
    AliTOFTenderSupply *tofTender = new AliTOFTenderSupply("TOFtender");
    tender->AddSupply(tofTender);
 }   

  //========= Attach TRD supply ======
  if (useTRD) {
    AliTRDTenderSupply *trdSupply=new AliTRDTenderSupply("TRDtender");

    trdSupply->SetLoadDeadChambersFromCDB();                    // Mask Bad chambers
    trdSupply->SetPIDmethod(AliTRDTenderSupply::k1DLQpid);
    trdSupply->SwitchOffGainCorrection();                       // Correction only on pass 1
    trdSupply->SetNormalizationFactor(0.12697,114737,130850);   // 1 otherwise
    trdSupply->SetRedoTRDMatching(kTRUE);
    tender->AddSupply(trdSupply);
  }  

  //========= Attach Primary Vertex supply ======
  if (useVTX) {
    tender->AddSupply(new AliVtxTenderSupply("PriVtxtender"));
  }  

  //========= Attach EMCAL supply ======
  if (useEmc) {
    AliEMCALTenderSupply *emcSupply = new AliEMCALTenderSupply("EmcTender");
    emcSupply->SetDefaults();
    tender->AddSupply(emcSupply);
  }  

  //========= Attach PID supply ======
  if (usePID) {
    AliPIDTenderSupply *pidSupply=new AliPIDTenderSupply("PIDtender");
    tender->AddSupply(pidSupply);
  }

  //================================================
  //              data containers
  //================================================

    //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("tender_event", AliESDEvent::Class(),
                           AliAnalysisManager::kExchangeContainer,"default_tender");
 
  //           connect containers
  mgr->ConnectInput  (tender,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (tender,  1, coutput1);
 
  return tender;
}
