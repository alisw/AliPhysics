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
  AliAnalysisTaskJetCorrel *task = new AliAnalysisTaskJetCorrel(Selector);
  mgr->AddTask(task);

  //create data containers
  // AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer(); // OLD WAY!
  AliAnalysisDataContainer *cinput = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");
  mgr->ConnectInput(task,0,cinput);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("JetCorrelHistos", TList::Class(),
			 AliAnalysisManager::kOutputContainer,"JetCorrelHistos.root");  
  mgr->ConnectOutput(task,0,coutput);
  
  return task;
}
