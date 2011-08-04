//=============================================//
const char* centralityEstimator = "V0M";
//=============================================//

//_________________________________________________________//
AliAnalysisTaskBF *AddTaskBalanceCentralityTrain(Double_t centrMin=0.,
						 Double_t centrMax=100.,
						 Double_t vertexZ=10.,
						 Double_t DCAxy=2.4,
						 Double_t DCAz=3.2,
						 Double_t ptMin=0.3,
						 Double_t ptMax=1.5,
						 Double_t etaMin=-0.8,
						 Double_t etaMax=0.8,
						 TString fileNameBase="AnalysisResults") {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);
  centralityName+="_";
  centralityName+=Form("vZ%.1f",vertexZ);
  centralityName+="_";
  centralityName+=Form("DCAxy%.1f",DCAxy);
  centralityName+="_";
  centralityName+=Form("DCAz%.1f",DCAz);
  centralityName+="_Pt";
  centralityName+=Form("%.1f",ptMin);
  centralityName+="-";
  centralityName+=Form("%.1f",ptMax);
  centralityName+="_Eta";
  centralityName+=Form("%.1f",etaMin);
  centralityName+="-";
  centralityName+=Form("%.1f",etaMax);

  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  gROOT->LoadMacro("./configBalanceFunctionAnalysis.C");
  AliBalance *bf  = 0;  // Balance Function object
  AliBalance *bfs = 0;  // shuffled Balance function object

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD");
    bfs = GetBalanceFunctionObject("ESD");
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD");
    bfs = GetBalanceFunctionObject("AOD");
  }
  else{
    bf  = GetBalanceFunctionObject("MC");
    bfs = GetBalanceFunctionObject("MC");
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBF *taskBF = new AliAnalysisTaskBF("TaskBF");
  taskBF->SetAnalysisObject(bf,bfs);
  taskBF->SetCentralityPercentileRange(centrMin,centrMax);
  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject();
    taskBF->SetAnalysisCutObject(trackCuts);
    taskBF->SetCentralityEstimator(centralityEstimator);
  }
  else if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);
  }

    // vertex cut (x,y,z)
    taskBF->SetVertexDiamond(.3,.3,vertexZ);

    // offline trigger selection (AliVEvent.h)
    // taskBF->UseOfflineTrigger(); // NOT usable 
    // with this only selected events are analyzed
    taskBF->SelectCollisionCandidates(AliVEvent::kMB);

    //bf->PrintAnalysisSettings();
    mgr->AddTask(taskBF);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2EbyE.outputBalanceFunctionAnalysis";
  AliAnalysisDataContainer *coutBalance = mgr->CreateContainer(Form("bfOutput_%s",centralityName.Data()), AliBalance::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(Form("listBF_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput(taskBF, 1, coutBalance);
  mgr->ConnectOutput(taskBF, 3, coutQA);
  mgr->ConnectOutput(taskBF, 4, coutBF);

  return taskBF;
}
