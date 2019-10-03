AliAnalysisTaskCTauPbPb* 
AddTaskCTau(Double_t min=0., Double_t max=90., 
TString name="cTau_0090", Bool_t isMC=kFALSE) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCTau", "No analysis manager to connect to.");
    return NULL;
  }  
  
  if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskCTau","This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskCTauPbPb *task = new AliAnalysisTaskCTauPbPb(name);
  task->SetCentrality(min,max);
  task->SetMC(isMC);
  mgr->AddTask(task);
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  if (isMC) name+="_mc";

  AliAnalysisDataContainer *coutput1 = 
     mgr->CreateContainer(name, TList::Class(),
     AliAnalysisManager::kOutputContainer, name+".root");
  mgr->ConnectOutput(task,1,coutput1);

  return task;
}   
