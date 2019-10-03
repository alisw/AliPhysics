//_________________________________________________________//
AliAnalysisTaskEfficiencyBF *AddTaskBalanceEfficiency(
						      TString  centralityEstimator="V0M",
						      Double_t centrMin=0.,
						      Double_t centrMax=90.,
						      Double_t vertexZ=10.,
						      TString fileNameBase="AnalysisResults"
						     ) {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");

  centralityName+=Form("%s_%.0f-%.0f_%.0f",centralityEstimator.Data(),centrMin,centrMax,vertexZ);
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");
  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTriggeredBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTriggeredBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";

 
  // Create the task, add it to manager and configure it.
  //===========================================================================


  AliAnalysisTaskEfficiencyBF *taskEfficiencyBF = new AliAnalysisTaskEfficiencyBF("TaskEfficiencyBF");

  // analysis mode
  taskEfficiencyBF->SetAnalysisMode("TPC");

  // centrality
  taskEfficiencyBF->SetCentralityEstimator(centralityEstimator);
  taskEfficiencyBF->SetCentralityPercentileRange(centrMin,centrMax);
  taskEfficiencyBF->SelectCollisionCandidates(AliVEvent::kMB);

  // vertex
  taskEfficiencyBF->SetVertexDiamond(.3,.3,vertexZ);

  // track cuts
  AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  cuts->SetRequireTPCStandAlone(kTRUE); // TPC only cuts!  
  taskEfficiencyBF->SetAnalysisCutObject(cuts);

  //analysis kinematic cuts
  taskEfficiencyBF->SetMinPt(0.3);
  taskEfficiencyBF->SetMaxPt(1.5);
  // taskEfficiencyBF->SetMinEta(-0.8);
  //taskEfficiencyBF->SetMaxEta(0.8);
  taskEfficiencyBF->SetEtaRange(-0.8,0.8,100,0.0,1.6, 64); //acceptance cuts
  taskEfficiencyBF->SetPtRange(0.1, 5.0, 49);  //acceptance cuts
  taskEfficiencyBF->SetPhiRange(0.0,360.,100, 180., 90);  //acceptance cuts
    
  // ADD the task
  //===========================================================================
  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskEfficiencyBF);

  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionEfficiencyAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutEfficiencyBF = mgr->CreateContainer(Form("listEfficiencyBF_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
 
  mgr->ConnectInput(taskEfficiencyBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEfficiencyBF, 1, coutQA);
  mgr->ConnectOutput(taskEfficiencyBF, 2, coutEfficiencyBF);

  return taskEfficiencyBF;
}
