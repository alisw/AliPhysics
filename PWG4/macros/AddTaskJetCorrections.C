AliAnalysisTaskJetCorrections * AddTaskJetCorrections()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetSpectrum", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetSpectrum", "This task requires an input event handler");
    return NULL;
   }
  
  AliAnalysisTaskJetCorrections * jetCorr = new  AliAnalysisTaskJetCorrections("Jet Corrections");
  
  jetCorr->SetBranchGen("jetsMC"); 
  jetCorr->SetBranchRec("jets");
  jetCorr->SetR(.5); 
  mgr->AddTask(jetCorr);
   
  AliAnalysisDataContainer *coutput1_Corr = mgr->CreateContainer("jetCorr", TList::Class(),AliAnalysisManager::kOutputContainer,"jetCorr.root");

   mgr->ConnectInput  (jetCorr, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (jetCorr, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (jetCorr,  1, coutput1_Corr );
   
   return jetCorr;
}


AliAnalysisTaskJetCorrections * AddTaskJetCorrections(AliAnalysisManager* mgr ,AliAnalysisDataContainer *cinput)
{
  if (!mgr) {
    ::Error("AddTaskJetSpectrum", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetSpectrum", "This task requires an input event handler");
    return NULL;
   }
  
  AliAnalysisTaskJetCorrections * jetCorr = new  AliAnalysisTaskJetCorrections("Jet Corrections");
  
  jetCorr->SetBranchGen("jetsMC"); 
  jetCorr->SetBranchRec("jets");
  jetCorr->SetR(.5); 
  mgr->AddTask(jetCorr);

  AliAnalysisDataContainer * coutpu0 = mgr->CreateContainer("coutpu0", TTree::Class(),
				  AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer *coutput1_jetCorr = mgr->CreateContainer("jetCorr", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_jetCorr",AliAnalysisManager::GetCommonFileName()));

   mgr->ConnectInput  (jetCorr, 0, cinput);
   mgr->ConnectOutput (jetCorr, 0, coutpu0);
   mgr->ConnectOutput (jetCorr,  1, coutput1_Corr );
   
   return jetCorr;
}
