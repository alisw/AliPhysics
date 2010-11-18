AliAnalysisTaskChargeFluctuations *AddTaskChargeFluctuations() {
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
    ::Error("AddTaskChargeFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  //gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/macros/configChargeFluctuationsAnalysis.C");
  gROOT->LoadMacro("configChargeFluctuationsAnalysis.C");

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskChargeFluctuations *taskChargeFluctuations = new AliAnalysisTaskChargeFluctuations("TaskChargeFluctuations");
  mgr->AddTask(taskChargeFluctuations);
  taskChargeFluctuations->SetAnalysisObject(bf);
  taskChargeFluctuations->SetVertexDiamond(0.3,0.3,10.);

  if (type=="ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject();
    taskChargeFluctuations->SetAnalysisCutObject(trackCuts);
  }

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2EbyE.outputChargeFluctuationsAnalysis.root";
  AliAnalysisDataContainer *coutChargeFluctuations = mgr->CreateContainer("outputList", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskChargeFluctuations, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskChargeFluctuations, 1, coutChargeFluctuations);

  // Return task pointer at the end
  return taskChargeFluctuations;
}
