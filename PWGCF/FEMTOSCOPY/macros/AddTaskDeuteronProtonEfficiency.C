
AliAnalysisTask *AddTaskDeuteronProtonEfficiency(){


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if(!mgr)
    {
      
      Error("AddTaskDeuteronProtonEfficiency.C","No analysis manager found.");
      return 0;

    }


  bool hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  if(!hasMC)
    {
        
      Error("AddTaskDeuteronProtonEfficiency.C","No MC Handler found.");
      //task->SelectCollisionCandidates(AliVEvent::kAny);
  //    return 0;

    }




  AliAnalysisTaskDeuteronProtonEfficiency *task = new AliAnalysisTaskDeuteronProtonEfficiency("DeuteronProtonEfficiency");


  mgr->AddTask(task);


  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput_dp_efficiency = mgr->CreateContainer("ProtonDeuteronEfficiency", TList::Class(), AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput_dp_efficiency);

  return task;



}
