//=============================================//
const char* centralityEstimator = "TRK";
//=============================================//

//_________________________________________________________//
AliEbyEFluctuationAnalysisTask *AddTaskFluctuations(const char* analysisType = "ESD", 
						    const char *analysisMode = "TPC") {
  // Creates a fluctuations analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliEbyEFluctuationAnalysisTask *taskFluctuations = new AliEbyEFluctuationAnalysisTask("FluctuationsTask");
  mgr->AddTask(taskFluctuations);
  taskFluctuations->SetAnalysisType(analysisType);
  if(analysisType == "ESD") {
    taskFluctuations->SetAnalysisMode(analysisMode);
    taskFluctuations->SetVertexDiamond(3.,3.,10.);
    taskFluctuations->SetCentralityEstimator(centralityEstimator);
    taskFluctuations->SetCentralityBins20();
    taskFluctuations->SelectCollisionCandidates(AliVEvent::kMB);
  }

  TString type = analysisType;
  if (analysisType=="ESD") {
    gROOT->LoadMacro("configFluctuationsAnalysis.C");
    AliESDtrackCuts *trackCuts = GetTrackCutsObject();
    taskFluctuations->SetAnalysisCutObject(trackCuts);
  }
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":outputFluctuationAnalysis.root";
  AliAnalysisDataContainer *coutFA = mgr->CreateContainer("fluctuationsOutput", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskFluctuations, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskFluctuations, 1, coutFA);

  // Return task pointer at the end
  return taskFluctuations;
}
