AliAnalysisTaskProtons *AddTaskBalanceFunction() {
  // Creates a proton analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskProtons", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBalanceFunction", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/macros/configBalanceFunctionAnalysis.C");
  AliBalance *bf = 0;
  if (type=="ESD") bf = GetBalanceFunctionObject("ESD");
  else if (type=="AOD") bf = GetBalanceFunctionObject("AOD");
  else bf = GetBalanceFunctionObject("MC");

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBF *taskBF = new AliAnalysisTaskBF("TaskBF");
  mgr->AddTask(taskBF);
  taskBF->SetAnalysisObject(bf);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2EbyE.outputBalanceFunctionAnalysis.root";
  AliAnalysisDataContainer *cout_bf = mgr->CreateContainer("bfOutput", AliBalance::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 0, cout_bf);

  // Return task pointer at the end
  return taskBF;
}
