AliAnalysisTaskJetCorrel *AddTaskJetCorrel(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    ::Error("AddTaskJetCorrel", "No analysis manager found.");
    return NULL;
  }
  if(!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetCorrel", "This task requires an input event handler.");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(inputDataType.CompareTo("ESD")!=0) {
    ::Error("AddTaskJetCorrel", "This task uses ESD as input.");
    return NULL;
  }

  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/ConfigJetCorrel.C");
  AliJetCorrelSelector* Selector = ConfigJetCorrel();
  AliAnalysisTaskJetCorrel *jetcorrel = new AliAnalysisTaskJetCorrel(Selector);
  mgr->AddTask(jetcorrel);
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  AliPhysicsSelection* physSele = physSelTask->GetPhysicsSelection();
  jetcorrel->SelectCollisionCandidates(); 

  //create data containers
  AliAnalysisDataContainer *cinput = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");
  mgr->ConnectInput(jetcorrel,0,cinput);
  AliAnalysisDataContainer *coutput_JetCorrel = mgr->CreateContainer("JetCorrelHistos", TList::Class(),
	AliAnalysisManager::kOutputContainer, Form("%s:JetCorrel",AliAnalysisManager::GetCommonFileName()));  
  mgr->ConnectOutput(jetcorrel,0,coutput_JetCorrel);
  
  return jetcorrel;
}
