AliAnalysisTaskCharmDecayTracks* AddTaskCharmDecayTracks(Int_t pdgCode=421,
							 Int_t optMeth=0,
							 Int_t filterMask = 16){
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCharmDecayTracks", "No analysis manager to connect to.");
  }
  
  AliAnalysisTaskCharmDecayTracks *dTask = new AliAnalysisTaskCharmDecayTracks();
  dTask->SetReadMC(kTRUE);
  dTask->SetFilterMask(filterMask);
  dTask->SetSelectedHadron(pdgCode);
  if(optMeth==0) dTask->SetUseCharmedHadronsFromKine();
  else dTask->SetUseCandidatesFromDeltaAOD();

  // Create containers for input/output
  TString baseName="CharmDecayTracks";
  TString containerStr="";
  if(pdgCode==421) containerStr="Dzero";
  else if(pdgCode==411) containerStr="Dplus";
  else if(pdgCode==431) containerStr="Ds";
  else if(pdgCode==4122) containerStr="Lc";
  TString inname = Form("cinput%s%s",baseName.Data(),containerStr.Data());
  TString outname = Form("coutput%s%s",baseName.Data(),containerStr.Data());
  TString treename = Form("coutput%sTree%s",baseName.Data(),containerStr.Data());
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),
                                                          AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":%s%s",baseName.Data(),containerStr.Data());
  
  
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputfile.Data());
  AliAnalysisDataContainer *coutputTree = mgr->CreateContainer(treename,TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               outputfile.Data());
  
  coutputTree->SetSpecialOutput();

  mgr->ConnectInput(dTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dTask,1,coutput);
  mgr->ConnectOutput(dTask,2,coutputTree);
  
  return dTask;

}
