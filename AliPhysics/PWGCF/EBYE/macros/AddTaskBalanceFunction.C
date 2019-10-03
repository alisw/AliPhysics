AliAnalysisTaskBF *AddTaskBalanceFunction() {
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
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/macros/configBalanceFunctionAnalysis.C");
  AliBalance *bf = 0;
  if (type=="ESD") bf = GetBalanceFunctionObject("ESD");
  else if (type=="AOD") bf = GetBalanceFunctionObject("AOD");
  else bf = GetBalanceFunctionObject("MC");

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBF *taskBF = new AliAnalysisTaskBF("TaskBF");
  mgr->AddTask(taskBF);
  taskBF->SetAnalysisObject(bf);
  taskBF->SetVertexDiamond(0.3,0.3,10.);

  if (type=="ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject();
    taskBF->SetAnalysisCutObject(trackCuts);
  }

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionAnalysis.root";
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer("bfOutput", AliBalance::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer("listQA", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutBF);
  mgr->ConnectOutput(taskBF, 2, coutQA);

  // Return task pointer at the end
  return taskBF;
}
