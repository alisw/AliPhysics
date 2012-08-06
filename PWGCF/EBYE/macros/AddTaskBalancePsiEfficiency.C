//_________________________________________________________//
AliAnalysisTaskEfficiencyBFPsi *AddTaskBalancePsiEfficiency(TString centralityEstimator="V0M", Double_t centrMin=0., Double_t centrMax=80., Double_t vertexZ=10., TString fileNameBase="AnalysisResults") {
  // Create a balance function analysis task and add it to the analysis manager
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");

  centralityName+=Form("%s_%.0f-%.0f_%.0f",centralityEstimator.Data(),centrMin,centrMax,vertexZ);
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBalancePsiEfficiency", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBalancePsiEfficiency", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";

 
  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskEfficiencyBFPsi *taskEfficiencyBFPsi = new AliAnalysisTaskEfficiencyBFPsi("TaskEfficiencyBFPsi");

  // analysis mode
  taskEfficiencyBFPsi->SetAnalysisMode("TPC");

  // centrality
  taskEfficiencyBFPsi->SetCentralityEstimator(centralityEstimator);
  taskEfficiencyBFPsi->SetCentralityPercentileRange(centrMin,centrMax);
  taskEfficiencyBFPsi->SelectCollisionCandidates(AliVEvent::kMB);

  // vertex
  taskEfficiencyBFPsi->SetVertexDiamond(.3,.3,vertexZ);

  // track cuts
  AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  cuts->SetRequireTPCStandAlone(kTRUE); // TPC only cuts!  
  taskEfficiencyBFPsi->SetAnalysisCutObject(cuts);

  //analysis kinematic cuts
  taskEfficiencyBFPsi->SetMinPt(0.2);
  taskEfficiencyBFPsi->SetMaxPt(20.0);
  taskEfficiencyBFPsi->SetMaxEta(0.8);
  taskEfficiencyBFPsi->SetEtaRangeMax(0.8); //acceptance cuts
  taskEfficiencyBFPsi->SetPtRangeMin(0.1);  //acceptance cuts
  taskEfficiencyBFPsi->SetPtRangeMax(5.0);  //acceptance cuts
  taskEfficiencyBFPsi->SetPhiRangeMin(0.);  //acceptance cuts
  taskEfficiencyBFPsi->SetPhiRangeMax(6.28);//acceptance cuts
  
  // ADD the task
  //===========================================================================
  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskEfficiencyBFPsi);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionPsiEfficiencyAnalysis";
  AliAnalysisDataContainer *coutQAPsi = mgr->CreateContainer(Form("listQAPsi_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutEfficiencyBFPsi = mgr->CreateContainer(Form("listEfficiencyBFPsi_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
 
  mgr->ConnectInput(taskEfficiencyBFPsi, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEfficiencyBFPsi, 1, coutQAPsi);
  mgr->ConnectOutput(taskEfficiencyBFPsi, 2, coutEfficiencyBFPsi);

  return taskEfficiencyBFPsi;
}
