AliAnalysisTask *AddTaskFindableHypertriton3(TString suffix = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFindableHypertriton3", "No analysis manager found.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFindableHypertriton3", "This task requires an input event handler");
    return 0x0;
  }

  // Create and configure the task
  TString tskname = "AliAnalysisTaskFindableHypertriton3";
  tskname.Append(suffix.Data());
  AliAnalysisTaskFindableHypertriton3 *task = new AliAnalysisTaskFindableHypertriton3(tskname.Data());
  mgr->AddTask(task);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  TString outputFileNameTree = outputFileName + ":FindableTree";

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *hypCont1 = mgr->CreateContainer("Hyp3FindTask_summary", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  AliAnalysisDataContainer *hypCont2 = mgr->CreateContainer("cTreeHyperTriton", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileNameTree);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, hypCont1);
  mgr->ConnectOutput(task, 2, hypCont2);

  return task;
}
