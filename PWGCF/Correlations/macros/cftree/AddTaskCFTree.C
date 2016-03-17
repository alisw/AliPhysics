AliAnalysisTaskCFTree *AddTaskCFTree(
    UInt_t selectionBit = AliVEvent::kMB,
    Double_t zVtxMax = 10,
    Double_t etaMax = 1,
    Double_t ptMin = 0.15,
    const char* outputFileName = 0,
    const char* folderName = "CorrelationTree")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhiCorrelations", "No analysis manager to connect to.");
    return NULL;
  }  
  
  AliAnalysisTaskCFTree* ana = new AliAnalysisTaskCFTree();
  ana->SetEventSelectionBit(selectionBit);
  ana->SetZVertex(zVtxMax);
  ana->SetTrackEtaCut(etaMax);
  ana->SetPtMin(ptMin);
  mgr->AddTask(ana);
  
  if (!outputFileName) outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histos", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1);
  mgr->ConnectOutput (ana, 2, coutput2);
  return ana;
}

