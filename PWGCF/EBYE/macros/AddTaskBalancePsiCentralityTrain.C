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
AliAnalysisTaskBFPsi *AddTaskBalancePsiCentralityTrain(Double_t centrMin=0.,
						       Double_t centrMax=100.,
						       Bool_t gRunShuffling=kFALSE,
						       Bool_t gRunMixing=kTRUE,
						       Bool_t gRunMixingWithEventPlane=kFALSE,
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
						       Bool_t bResonancesCut = kTRUE,
						       Bool_t bHBTcut = kTRUE,
						       Bool_t bConversionCut = kTRUE,
						       Bool_t bMomentumDifferenceCut = kTRUE,
						       Int_t AODfilterBit = 128,
						       Bool_t bCentralTrigger = kFALSE,
						       TString fileNameBase="AnalysisResults",
						       TString fArgEventClass="Centrality",
						       TString analysisTypeUser="AOD",
						       Bool_t bVertexBinning=kTRUE) {
  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
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
  if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
  
  // to set the analysis type manually
  if(analysisTypeUser != ""){
    analysisType = analysisTypeUser;
    ::Info("AddTaskBF",Form("Analysis Type manually set to %s",analysisType.Data()));
  }

  // for local changed BF configuration
  //gROOT->LoadMacro("./configBalanceFunctionPsiAnalysis.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGCF/EBYE/macros/configBalanceFunctionPsiAnalysis.C");
  AliBalancePsi *bf  = 0;  // Balance Function object
  AliBalancePsi *bfs = 0;  // shuffled Balance function object
  AliBalancePsi *bfm = 0;  // mixing Balance function object

  //maximum Delta eta range
  Double_t deltaEtaMax=TMath::Abs(etaMax-etaMin);

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("ESD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("AOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
  }
  else if (analysisType=="MC"){
    bf  = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("MC",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
  }
  else if (analysisType=="MCAOD"){
    bf  = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunShuffling) bfs = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kTRUE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
    if(gRunMixing)    bfm = GetBalanceFunctionObject("MCAOD",centralityEstimator,centrMin,centrMax,kFALSE,bResonancesCut,bHBTcut,bConversionCut,bMomentumDifferenceCut,fArgEventClass,deltaEtaMax,bVertexBinning);
  }
  else{
    ::Error("AddTaskBF", "analysis type NOT known.");
    return NULL;
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBFPsi *taskBF = new AliAnalysisTaskBFPsi("TaskBFPsi");
  
  //Event characteristics scheme
  taskBF->SetEventClass(fArgEventClass);
 
  //++++++++++++++++++++++
  TString correctionFileName = "$ALICE_ROOT/PWGCF/EBYE/BalanceFunctions/Corrections/CorrectionMaps.root"; //to put the path for the correction maps
  TFile *fCorrectionMatrix  = TFile::Open(correctionFileName.Data());
  if(!fCorrectionMatrix){
    Printf("WARNING CORRECTION histogram file not found");
    return NULL;
  }
  
  taskBF->SetInputCorrection(correctionFileName.Data(),"");

  //+++++++++++++++++++++

  taskBF->SetAnalysisObject(bf);
  if(gRunShuffling) taskBF->SetShufflingObject(bfs);
  if(gRunMixing){
    taskBF->SetMixingObject(bfm);
    taskBF->SetMixingTracks(50000);
    if(gRunMixingWithEventPlane){
      taskBF->SetMixingWithEventPlane(gRunMixingWithEventPlane);
    }
  }

  taskBF->SetCentralityPercentileRange(centrMin,centrMax);
  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject(ptMin,ptMax,etaMin,etaMax,maxTPCchi2,DCAxy,DCAz,minNClustersTPC);
    taskBF->SetAnalysisCutObject(trackCuts);
    if(kUsePID) {
      if(kUseBayesianPID)
	taskBF->SetUseBayesianPID(gMinAcceptedProbability);
      else if(kUseNSigmaPID)
	taskBF->SetUseNSigmaPID(nSigmaMax);
      taskBF->SetParticleOfInterest(AliAnalysistaskBFPsi::kProton);
      taskBF->SetDetectorUsedForPID(AliAnalysisTaskBFPsi::kTOFpid);
    }
  }
  else if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);

    // set extra DCA cuts (-1 no extra cut)
    taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);

    // set extra TPC chi2 / nr of clusters cut
    taskBF->SetExtraTPCCutsAOD(maxTPCchi2, minNClustersTPC);
    
  }
  else if(analysisType == "MC") {
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax); 
    taskBF->SetImpactParameterRange(centrMin,centrMax);
  }
  else if(analysisType == "MCAOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(AODfilterBit);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);    
  }

  // offline trigger selection (AliVEvent.h)
  // taskBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
  // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
  // documentation in https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PWG1EvSelDocumentation
  if(bCentralTrigger) taskBF->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  else                taskBF->SelectCollisionCandidates(AliVEvent::kMB);

  // centrality estimator (default = V0M)
  taskBF->SetCentralityEstimator(centralityEstimator);
  
  // vertex cut (x,y,z)
  taskBF->SetVertexDiamond(3.,3.,vertexZ);
  


  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskBF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCFEbyE.outputBalanceFunctionPsiAnalysis";
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(Form("listQAPsi_%.0f-%.0f_Bit%d_%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(Form("listBFPsi_%.0f-%.0f_Bit%d_%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunShuffling) AliAnalysisDataContainer *coutBFS = mgr->CreateContainer(Form("listBFPsiShuffled_%.0f-%.0f_Bit%d_%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(gRunMixing) AliAnalysisDataContainer *coutBFM = mgr->CreateContainer(Form("listBFPsiMixed_%.0f-%.0f_Bit%d_%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  if(kUsePID) AliAnalysisDataContainer *coutQAPID = mgr->CreateContainer(Form("listQAPIDPsi_%.0f-%.0f_Bit%d_%s",centrMin,centrMax,AODfilterBit,centralityEstimator.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutQA);
  mgr->ConnectOutput(taskBF, 2, coutBF);
  if(gRunShuffling) mgr->ConnectOutput(taskBF, 3, coutBFS);
  if(gRunMixing) mgr->ConnectOutput(taskBF, 4, coutBFM);
  if(kUsePID && analysisType == "ESD") mgr->ConnectOutput(taskBF, 5, coutQAPID);

  return taskBF;
}
