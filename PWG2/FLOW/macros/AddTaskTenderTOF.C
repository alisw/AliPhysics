AliAnalysisTask *AddTaskTenderTOF(Float_t tofres = 80,Bool_t corrExpTimes=kFALSE){
  //get the current analysis manager
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
  tender->SetCheckEventSelection(kFALSE);
  //tender->SetDefaultCDBStorage("raw://");
  tender->SetDefaultCDBStorage("alien://folder=/alice/data/2010/OCDB");
  mgr->AddTask(tender);
  
  //========= Attach TOF supply ======
  AliTOFTenderSupply *TOFtender = new AliTOFTenderSupply("TOFtender");
  TOFtender->SetTOFres(tofres);
  TOFtender->SetCorrectExpTimes(corrExpTimes);
  //TOFtender->SetTheorExpTimes(kTRUE);
  tender->AddSupply(TOFtender);
  
    //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("tender_event", AliESDEvent::Class(),
                           AliAnalysisManager::kExchangeContainer,"default_tender");
 
  //           connect containers
  mgr->ConnectInput  (tender,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (tender,  1, coutput1);
 
  return tender;
}
