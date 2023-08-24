AliAnalysisSigmaBarCharged* AddTaskPHOSSigmaBar(const char* name = "PHOSSigmaBar", const bool isMC = false,
                                                Float_t eMinCut = 0.8, Float_t tCut = 25.e-9, int trackbits=21)
{
  // Add a task AliAnalysisSigmaBarCharged to the analysis train
  // Author: Pavel Gordeev, D.Peresunko, NRC "Kurchatov institute"

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSTagging", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSTagging", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisSigmaBarCharged* task = new AliAnalysisSigmaBarCharged(name);
  task->SetMC(isMC);
  task->SetMinNbarEnergy(eMinCut);
  task->SetClusterTOF(tCut);
  task->SetTrackBits(trackbits) ;
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  AliAnalysisDataContainer* coutput1 =
    mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, "histos.root");
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}
