// now in options
//=============================================//
//const char* centralityEstimator = "V0M";
//const char* centralityEstimator = "CL1";
//const char* centralityEstimator = "TRK";
//=============================================//
//Bool_t gRunShuffling = kFALSE;
//Bool_t gRunShuffling = kTRUE;
//=============================================//

//PID config
Bool_t kUseNSigmaPID = kFALSE;
Double_t nSigmaMax = 3.0;
Bool_t kUseBayesianPID = kTRUE;
Double_t gMinAcceptedProbability = 0.7;

//_________________________________________________________//
AliAnalysisTaskEventMixingBF *AddTaskBalanceEventMixing(Double_t centrMin=0.,
						 Double_t centrMax=100.,
						 Bool_t gRunShuffling=kFALSE,
						 TString centralityEstimator="V0M",
						 Double_t vertexZ=10.,
						 Double_t DCAxy=-1,
						 Double_t DCAz=-1,
						 Double_t ptMin=0.3,
						 Double_t ptMax=1.5,
						 Double_t etaMin=-0.8,
						 Double_t etaMax=0.8,
						 Double_t maxTPCchi2 = -1, 
						 Int_t minNClustersTPC = -1,
						 Bool_t kUsePID = kFALSE,
						 Int_t AODfilterBit = 128,
						 Bool_t bCentralTrigger = kFALSE,
						 TString fileNameBase="AnalysisResults") {

  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString centralityName("");
  centralityName+=Form("%.0f",centrMin);
  centralityName+="-";
  centralityName+=Form("%.0f",centrMax);
  centralityName+="_";
  centralityName+=Form("%s",centralityEstimator.Data());
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
  centralityName+="_Chi";
  centralityName+=Form("%.1f",maxTPCchi2);
  centralityName+="_nClus";
  centralityName+=Form("%d",minNClustersTPC);
  centralityName+="_Bit";
  centralityName+=Form("%d",AODfilterBit);
  if(bCentralTrigger)   centralityName+="_withCentralTrigger";





  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEventMixingBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventMixingBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";

  // for local changed EventMixingBF configuration
  //gROOT->LoadMacro("./configBalanceFunctionAnalysisEventMixing.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGCF/EBYE/macros/configBalanceFunctionAnalysisEventMixing.C");
  AliBalanceEventMixing *bf  = 0;  // Balance Function object
  AliBalanceEventMixing *bfs = 0;  // shuffled Balance function object

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kTRUE);
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kTRUE);
  }
  else if (analysisType=="MC"){
    bf  = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kTRUE);
  }
  else{
    ::Error("AddTaskEventMixingBF", "analysis type NOT known.");
    return NULL;
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskEventMixingBF *taskEventMixingBF = new AliAnalysisTaskEventMixingBF("TaskEventMixingBF");
  taskEventMixingBF->SetAnalysisObject(bf);
  if(gRunShuffling) taskEventMixingBF->SetShufflingObject(bfs);

  taskEventMixingBF->SetCentralityPercentileRange(centrMin,centrMax);
  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject(ptMin,ptMax,etaMin,etaMax,maxTPCchi2,DCAxy,DCAz,minNClustersTPC);
    taskEventMixingBF->SetAnalysisCutObject(trackCuts);
    if(kUsePID) {
      if(kUseBayesianPID)
	taskEventMixingBF->SetUseBayesianPID(gMinAcceptedProbability);
      else if(kUseNSigmaPID)
	taskEventMixingBF->SetUseNSigmaPID(nSigmaMax);
      taskEventMixingBF->SetParticleOfInterest(AliAnalysistaskEventMixingBF::kProton);
      taskEventMixingBF->SetDetectorUsedForPID(AliAnalysisTaskEventMixingBF::kTOFpid);
    }
  }
  else if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskEventMixingBF->SetAODtrackCutBit(AODfilterBit);
    taskEventMixingBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);

    // set extra DCA cuts (-1 no extra cut)
    taskEventMixingBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskEventMixingBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);
    
  }
  else if(analysisType == "MC") {
    taskEventMixingBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax); 
  }

  // offline trigger selection (AliVEvent.h)
  // taskEventMixingBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
  // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
  // documentation in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWG1EvSelDocumentation
  if(bCentralTrigger) taskEventMixingBF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else                taskEventMixingBF->SelectCollisionCandidates(AliVEvent::kMB);

  // centrality estimator (default = V0M)
  taskEventMixingBF->SetCentralityEstimator(centralityEstimator);
  
  // vertex cut (x,y,z)
  taskEventMixingBF->SetVertexDiamond(.3,.3,vertexZ);
  


  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskEventMixingBF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQA_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutEventMixingBF = mgr->CreateContainer(Form("listEventMixingBF_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunShuffling) AliAnalysisDataContainer *coutEventMixingBFS = mgr->CreateContainer(Form("listEventMixingBFShuffled_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(kUsePID) AliAnalysisDataContainer *coutQAPID = mgr->CreateContainer(Form("listQAPID_%s",centralityName.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskEventMixingBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEventMixingBF, 1, coutQA);
  mgr->ConnectOutput(taskEventMixingBF, 2, coutEventMixingBF);
  if(gRunShuffling) mgr->ConnectOutput(taskEventMixingBF, 3, coutEventMixingBFS);
  if(kUsePID && analysisType == "ESD") mgr->ConnectOutput(taskEventMixingBF, 4, coutQAPID);

  return taskEventMixingBF;
}
