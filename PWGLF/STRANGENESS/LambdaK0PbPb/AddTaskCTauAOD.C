AliAnalysisTaskCTauPbPbaod* 
AddTaskCTauAOD(Double_t min=0., Double_t max=90., Double_t cpa=0.9975, 
Double_t dca=1.5, Double_t cr=70., Double_t crfd=0.8, Double_t pv=0.1, 
TString name="cTau_0090aod", Bool_t isMC=kFALSE) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskCTauAOD", "No analysis manager to connect to.");
    return NULL;
  }  
  
  if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskCTauAOD","This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskCTauPbPbaod *task = new AliAnalysisTaskCTauPbPbaod(name);
  task->SetCentrality(min,max);
  task->SetMC(isMC);

  task->SetCosPA(cpa);
  task->SetDtrDCA(dca);
  task->SetTPCrows(cr);
  task->SetTPCratio(crfd);
  task->SetPrimDCA(pv);

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
