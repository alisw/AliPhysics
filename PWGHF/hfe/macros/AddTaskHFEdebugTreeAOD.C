AliAnalysisTask *AddTaskHFEdebugTreeAOD(){

  // libraries in case
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr->GetInputEventHandler()) {
    printf("AddTask_hfe_HFEdebugTreeTaskAOD", "This task requires an input event handler");
    return NULL;
  }

  AliHFEdebugTreeTaskAOD *task = new AliHFEdebugTreeTaskAOD("HFEdebugTreeCreator");
  task->SetFileName("HFEdebug.root");
  task->SetMinNclustersTPC(30);
  task->SetMinNclustersITS(2);
  task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);

   
  mgr->AddTask(task);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += "debugtreeaod";

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(task,1, mgr->CreateContainer("debugtreeaod", TTree::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(task,0, cinput );    

  return task;
  
}
