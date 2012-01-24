//=============================================//
const char* centralityEstimator = "TRK";
//=============================================//

//_________________________________________________________//
void AddTaskFluctuationsCentralityTrain(const char* analysisType = "ESD", 
					const char *analysisMode = "TPC",
					Float_t centrMin=0.,
					Float_t centrMax=100.,
					TString fileNameBase="AnalysisResults") {
  // Creates a fluctuations analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);

  TString fileName(fileNameBase);
  fileName.Append(".root");

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
  AliEbyEFluctuationAnalysisTaskTrain *taskFluctuations = new AliEbyEFluctuationAnalysisTaskTrain("FluctuationsTask");
  taskFluctuations->SetCentralityPercentileRange(centrMin,centrMax);
  taskFluctuations->SetAnalysisType(analysisType);
  if(analysisType == "ESD") {
    taskFluctuations->SetAnalysisMode(analysisMode);
    taskFluctuations->SetVertexDiamond(3.,3.,10.);
    taskFluctuations->SetCentralityEstimator(centralityEstimator);
    //taskFluctuations->SetCentralityBins20();
    taskFluctuations->SelectCollisionCandidates(AliVEvent::kMB);
  }

  TString type = analysisType;
  //if (analysisType=="ESD") {
  //gROOT->LoadMacro("configFluctuationsAnalysis.C");
  //AliESDtrackCuts *trackCuts = GetTrackCutsObject();
  //taskFluctuations->SetAnalysisCutObject(trackCuts);
  //}
  
  mgr->AddTask(taskFluctuations);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutputFA = mgr->CreateContainer(Form("cobjFA_%s",centralityName.Data()), 
							     TList::Class(),AliAnalysisManager::kOutputContainer,fileName.Data()); 
  mgr->ConnectInput(taskFluctuations,0, mgr->GetCommonInputContainer()); 
  mgr->ConnectOutput(taskFluctuations,1,coutputFA); 

  return;
}
