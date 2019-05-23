AliAnalysisTaskNtuplizer* AddTaskNtuplizer()
{
   // Creates and configures AliAnalysisTaskNtuplizer class instance.

   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskNtuplizer", "mgr is NULL, task not added");
      return NULL;
   }

    if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPIDCorr", "This task requires an input event handler");
    return NULL;
  }

   TString type = mgr->GetInputEventHandler()->GetDataType();
   

   AliAnalysisTaskNtuplizer* task = new AliAnalysisTaskNtuplizer("TaskCheckData");
   // task->SelectCollisionCandidates(trg);
   mgr->AddTask(task);

   TString outpath = AliAnalysisManager::GetCommonFileName();
   outpath += ":Container";

   // output containers
   AliAnalysisDataContainer* out1 = mgr->CreateContainer("histos",
      TObjArray::Class(), AliAnalysisManager::kOutputContainer, outpath);

   AliAnalysisDataContainer* out2 = mgr->CreateContainer("TTree",
      TTree::Class(), AliAnalysisManager::kOutputContainer, outpath);

   out2->SetSpecialOutput();  // file-resident

   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, out1);
   mgr->ConnectOutput(task, 2, out2);

   return task;
}
