//______________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskReducedTreeNuclei (TString prod="LHC13b")  {

   //Get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("No analysis manager found.");
      return 0;
   }

   Int_t collSystem = 0;
   if(!prod.CompareTo("LHC10b") || !prod.CompareTo("LHC10c") || !prod.CompareTo("LHC10d") ||
      !prod.CompareTo("LHC10e") || !prod.CompareTo("LHC10f")) collSystem = 1;   // minimum bias pp
   if(!prod.CompareTo("LHC13b") || !prod.CompareTo("LHC13c")|| !prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")) collSystem = 2;   // minimum bias p-Pb
   if(!prod.CompareTo("LHC16l")) collSystem = 3;    // minimum bias pp and Pb-Pb

   AliAnalysisTaskReducedTreeNuclei *task = new AliAnalysisTaskReducedTreeNuclei("TaskNuclei");
   if(collSystem==1) task->SelectCollisionCandidates(AliVEvent::kMB);
   if(collSystem==2) task->SelectCollisionCandidates(AliVEvent::kINT7);
   if(collSystem==3) task->SelectCollisionCandidates(AliVEvent::kMB+AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
   mgr->AddTask(task);

   TString Filename = "ReducedTreeNuclei.root"; // mgr->GetCommonFileName();
   AliAnalysisDataContainer *cOutput = mgr->CreateContainer("ReducedTreeNuclei", TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   AliAnalysisDataContainer *cQA = mgr->CreateContainer("QAHistograms", TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
   mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, cOutput);
   mgr->ConnectOutput(task, 2, cQA);

   return task;
}
//______________________________________________________________________________________________________________________________________________________________________________________






