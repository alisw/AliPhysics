void AddTask_GammaPythia(Double_t maxpT = 20,
			 Int_t doMultStudies = 1)
{

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_GammaPythia", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskGammaPythia *task=NULL;
  task= new AliAnalysisTaskGammaPythia("GammaPythia");
  task->SetMaxPt(maxpT);
  task->SetDoMultStudies(doMultStudies);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer("GammaPythia", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaPythia",AliAnalysisManager::GetCommonFileName()));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
