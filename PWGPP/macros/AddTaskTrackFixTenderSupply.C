//______________________________________________________________________________
AliAnalysisTask* AddTaskTrackFixTenderSupply(const char* passName,
                                             const char* objOADBpath="$OADB/PWGPP/data/CorrPTInv.root",
                                             const char* ocdb="raw://" )
{
  //adds a tender to fix the momenta of tracks.
  //passName has to be provided as the correction are reconstruction pass dependent
  //the corresponding containers in OADB are named after the pass "passName"
  //TODO: make the pass detection automatic
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  // Tender and supplies. Needs to be called for every event.
  AliTender *tender=new AliTender("AnalysisTender");
  tender->SetDefaultCDBStorage(ocdb);
  
  AliTrackFixTenderSupply* ptinvCor = new AliTrackFixTenderSupply("ptinvCorrSupply");
  //ptinvCor->SetDebugLevel(0);
  ptinvCor->SetOADBObjPath(objOADBpath);
  ptinvCor->SetOADBObjName(passName);
  
  tender->AddSupply(ptinvCor);
  mgr->AddTask(tender);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      "trackFixCorrectionTender", 
      AliESDEvent::Class(),
      AliAnalysisManager::kExchangeContainer,
      "default_tender");

  //           connect containers
  mgr->ConnectInput  (tender,  0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (tender,  1, coutput1);

  return tender;
}
