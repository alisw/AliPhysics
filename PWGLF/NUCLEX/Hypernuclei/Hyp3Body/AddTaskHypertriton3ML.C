AliAnalysisTask *AddTaskHypertriton3ML(Bool_t readMC = kFALSE, TString tskname = "Hypertriton3body", TString suffix = "") {

  // Creates, configures and attaches to the train the task hyper-triton 3 body 
  // decay study
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  ::Info("AddTaskHypertriton3ML", "Adding a new task with this settings readMC = %i", readMC);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHypertriton3ML", "No analysis manager to connect to.");
    return 0x0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHypertriton3ML", "This task requires an input event handler");
    return 0x0;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type.Contains("AOD")) {
    ::Error("AddTaskHypertriton3ML", "This task requires to run on ESD");
    return 0x0;
  }

  // Create and configure the task
  tskname.Append(suffix.Data());
  AliAnalysisTaskHypertriton3New *taskhyp = new AliAnalysisTaskHypertriton3New(true, tskname.Data());
  
  mgr->AddTask(taskhyp);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *hypCont1 = mgr->CreateContainer(Form("%s_summary", tskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  AliAnalysisDataContainer *hypCont2 = mgr->CreateContainer(Form("Hyp3Tree%s", suffix.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("Hyp3Tree.root:%s",suffix.Data()));

  mgr->ConnectInput(taskhyp, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskhyp, 1, hypCont1);
  mgr->ConnectOutput(taskhyp, 2, hypCont2);

  return taskhyp;
}
