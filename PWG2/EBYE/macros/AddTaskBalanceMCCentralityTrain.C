// now in options
//=============================================//
//const char* centralityEstimator = "V0M";
//const char* centralityEstimator = "CL1";
//const char* centralityEstimator = "TRK";
//=============================================//
//Bool_t gRunShuffling = kFALSE;
//Bool_t gRunShuffling = kTRUE;
//=============================================//
//_________________________________________________________//
AliAnalysisTaskBF *AddTaskBalanceMCCentralityTrain(Double_t centrMin=0.,
						   Double_t centrMax=100.,
						   Double_t impactParameterMin=0.,
						   Double_t impactParameterMax=20.,
						   Bool_t gRunShuffling=kFALSE,
						   Double_t vertexZ=10.,
						   Double_t ptMin=0.3,
						   Double_t ptMax=1.5,
						   Double_t etaMin=-0.8,
						   Double_t etaMax=0.8,
						   TF1* gAcceptanceParameterization = 0x0,
						   Int_t gPdgCode = -1,
						   TString fileNameBase="AnalysisResults") {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);

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
  TString analysisType = "MC";

  // for local changed BF configuration
  gROOT->LoadMacro("./configBalanceFunctionAnalysis.C");
  //gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/macros/configBalanceFunctionAnalysis.C");
  AliBalance *bf  = 0;  // Balance Function object
  AliBalance *bfs = 0;  // shuffled Balance function object

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD",centralityName.Data());
    if(gRunShuffling) bfs = GetBalanceFunctionObject("ESD",centralityName.Data(),kTRUE);
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD",centralityName.Data());
    if(gRunShuffling) bfs = GetBalanceFunctionObject("AOD",centralityName.Data(),kTRUE);
  }
  else if (analysisType=="MC"){
    bf  = GetBalanceFunctionObject("MC",centralityName.Data());
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MC",centralityName.Data(),kTRUE);
  }
  else{
    ::Error("AddTaskBF", "analysis type NOT known.");
    return NULL;
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBF *taskBF = new AliAnalysisTaskBF("TaskBF");
  taskBF->SetAnalysisObject(bf);
  if(gRunShuffling) taskBF->SetShufflingObject(bfs);

  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject(ptMin,ptMax,etaMin,etaMax,maxTPCchi2,DCAxy,DCAz,minNClustersTPC);
    taskBF->SetAnalysisCutObject(trackCuts);
    // centrality estimator (default = V0M)
    taskBF->SetCentralityEstimator(centralityEstimator);
    taskBF->SetCentralityPercentileRange(impactParameterMin,
					 impactParameterMax);
  }
  else if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(128);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    taskBF->SetCentralityEstimator(centralityEstimator);
    taskBF->SetCentralityPercentileRange(impactParameterMin,
					 impactParameterMax);

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);
    taskBF->SetCentralityEstimator(centralityEstimator);    
  }
  else if(analysisType == "MC") {
    Printf("********************ANALYSIS TYPE MC********************************");
    if(gAcceptanceParameterization)
      taskBF->SetAcceptanceParameterization(gAcceptanceParameterization);
    if(gPdgCode != -1)
      taskBF->SetPDGCode(gPdgCode);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax); 
    taskBF->SetImpactParameterRange(impactParameterMin,
				    impactParameterMax);
  }

  // offline trigger selection (AliVEvent.h)
  // taskBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
  // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
  // documentation in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWG1EvSelDocumentation
  //taskBF->SelectCollisionCandidates(AliVEvent::kMB);
    
  // vertex cut (x,y,z)
  taskBF->SetVertexDiamond(.3,.3,vertexZ);
  


  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskBF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2EbyE.outputBalanceFunctionAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(Form("listBF_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunShuffling) AliAnalysisDataContainer *coutBFS= mgr->CreateContainer(Form("listBFShuffled_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutQA);
  mgr->ConnectOutput(taskBF, 2, coutBF);
  if(gRunShuffling) mgr->ConnectOutput(taskBF, 3, coutBFS);

  return taskBF;
}
