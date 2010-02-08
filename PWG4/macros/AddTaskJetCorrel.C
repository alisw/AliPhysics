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

  using namespace JetCorrelHD;
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/ConfigJetCorrel.C");
  AliJetCorrelSelector* Selector = ConfigJetCorrel();
  AliAnalysisTaskJetCorrel *task = new AliAnalysisTaskJetCorrel(Selector);

  //add the task to the current analysis manager
  mgr->AddTask(task);
  //----------------------
  //create data containers
  //----------------------
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("JetCorrelHistos", TList::Class(),
			 AliAnalysisManager::kOutputContainer,"JetCorrelHistos.root");
  AliAnalysisDataContainer *dummy = mgr->CreateContainer("cdummy", TTree::Class(), 
							 AliAnalysisManager::kExchangeContainer);  
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,0,dummy);
  mgr->ConnectOutput(task,1,output);
  
  return task;
}
