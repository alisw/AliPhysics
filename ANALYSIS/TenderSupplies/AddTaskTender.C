AliAnalysisTask *AddTaskTender(Bool_t useV0=kFALSE){
  //get the current analysis manager
  Bool_t checkEvtSelection = useV0;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_tender_Tender", "No analysis manager found.");
    return 0;
  }
  // currently don't accept AOD input
  if (!mgr->GetInputEventHandler()->InheritsFrom(AliESDInputHandler::Class())) {
    Error("AddTask_tender_Tender","The analysis tender only works with ESD input!");
    return 0;
  }

  
  //========= Add tender to the ANALYSIS manager and set default storage =====
  AliTender *tender=new AliTender("AnalysisTender");
  tender->SetCheckEventSelection(checkEvtSelection);
  tender->SetDefaultCDBStorage("raw://");
  mgr->AddTask(tender);
  if (checkEvtSelection) {
     if (mgr->GetTasks()->First() != (TObject*)tender) {
        ::Error("When setting the tender to check the event selection, it has to be the first wagon ! Aborting.");
        return NULL;
     }
  }   
  
  //========= Attach VZERO supply ======
  if (useV0) {
     AliVZEROTenderSupply *vzeroSupply=new AliVZEROTenderSupply("VZEROtender");
     vzeroSupply->SetDebug(kFALSE);
     tender->AddSupply(vzeroSupply);
  }   
  //========= Attach TPC supply ======
  AliTPCTenderSupply *tpcSupply=new AliTPCTenderSupply("TPCtender");
  tpcSupply->SetDebugLevel(2);
  //tpcSupply->SetMip(50.);
  tender->AddSupply(tpcSupply);

  //========= Attach TOF supply ======
  AliTOFTenderSupply *tofTender = new AliTOFTenderSupply("TOFtender");
  tender->AddSupply(tofTender);
  
  //========= Attach TRD supply ======
  AliTRDTenderSupply *trdSupply=new AliTRDTenderSupply("TRDtender");
  tender->AddSupply(trdSupply);

  //========= Attach PID supply ======
  tender->AddSupply(new AliPIDTenderSupply("PIDtender"));

  //========= Attach Primary Vertex supply ======
  tender->AddSupply(new AliVtxTenderSupply("PriVtxtender"));
  
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
