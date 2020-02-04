
AliAnalysisTask* AddTask(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun, UInt_t kTriggerInt, Float_t minc, Float_t maxc, Bool_t CentFrameworkAliCen)
{

  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  // access methodh
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHighPtDeDx", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskHighPtDeDx", "This task requires an input event handler");
    return NULL;
  }  
  
  //
  // Add track filters, with Golden Cuts
  //
 
  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");

  AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  esdTrackCutsGolden->SetMinNCrossedRowsTPC(120);
  esdTrackCutsGolden->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsGolden->SetMaxChi2PerClusterITS(36);
  esdTrackCutsGolden->SetMaxFractionSharedTPCClusters(0.4);
  esdTrackCutsGolden->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCutsGolden->SetMaxDCAToVertexXY(3.0);
  
  trackFilterGolden->AddCuts(esdTrackCutsGolden);
  
  //old cuts without golden cut
  AliAnalysisFilter* trackFilter0 = new AliAnalysisFilter("trackFilter");
  Bool_t clusterCut = 0;
  Bool_t selPrimaries = kTRUE;
  AliESDtrackCuts* esdTrackCutsL0 = new AliESDtrackCuts;
  // TPC  
  if(clusterCut == 0)  esdTrackCutsL0->SetMinNClustersTPC(70);
  else if (clusterCut == 1) {
    esdTrackCutsL0->SetMinNCrossedRowsTPC(70);
    esdTrackCutsL0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  }
  else {
    AliWarningClass(Form("Wrong value of the clusterCut parameter (%d), using cut on Nclusters",clusterCut));
    esdTrackCutsL0->SetMinNClustersTPC(70);
  }

  esdTrackCutsL0->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsL0->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsL0->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCutsL0->SetRequireITSRefit(kTRUE);
  esdTrackCutsL0->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  if(selPrimaries) {
    esdTrackCutsL0->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  }
  esdTrackCutsL0->SetMaxDCAToVertexZ(2);
  esdTrackCutsL0->SetDCAToVertex2D(kFALSE);
  esdTrackCutsL0->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsL0->SetMaxChi2PerClusterITS(1e10);
  esdTrackCutsL0->SetMaxChi2TPCConstrainedGlobal(1e10);

  trackFilter0->AddCuts(esdTrackCutsL0);
  
  AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
  AliESDtrackCuts* esdTrackCutsTPC = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackFilterTPC->AddCuts(esdTrackCutsTPC);


  // Create the task and configure it 
  //========================================================================
  if(typerun==2){//pbpb: heavy ion analysis
    
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx;
    
    taskHighPtDeDx = 0;
    Char_t TaskName[256] = {0};
    sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc,maxc);
    taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx(TaskName);

    // TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    // taskHighPtDeDx->SetAnalysisType(type);
    // Run AOD even when filtered from LF_PbPb or LF_PbPb_MC

    //Set analysis details
    taskHighPtDeDx->SetAnalysisType("AOD");
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE); //for pbpb and ppb
    taskHighPtDeDx->SetAnalysisRun2(kFALSE); 
    taskHighPtDeDx->SetTrigger1(kTriggerInt); 
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetDebugLevel(0);
    
    //Set event details
    taskHighPtDeDx->SetCentFrameworkAliCen(CentFrameworkAliCen); //kTRUE:AliCentrality, kFALSE:AliMultSelection
    taskHighPtDeDx->SetCentDetector("V0M");
    taskHighPtDeDx->SetMinCent(minc);
    taskHighPtDeDx->SetMaxCent(maxc);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetContributorsVtxCut(0);
    taskHighPtDeDx->SetContributorsVtxSPDCut(0);
    taskHighPtDeDx->SetPileupCut(999.); //Correlation between global Zvtx and SPD Zvtx: use 999. for no cut
    taskHighPtDeDx->SetVtxR2Cut(10.); //use 10. for no cut
    
    //Set trigger particle tracks, v0s and daughter track details
    taskHighPtDeDx->SetMinPt(5.0); //trigger tracks
    taskHighPtDeDx->SetLowPtFraction(0.0); //keep x.x% of tracks below min pt
    taskHighPtDeDx->SetEtaCut(0.8); 
    taskHighPtDeDx->SetCrossedRowsCut(70.); //use 0. for no cut
    taskHighPtDeDx->SetCrossedOverFindableCut(0.8); //use 0. for no cut
    taskHighPtDeDx->SetCosPACut(0.97);
    taskHighPtDeDx->SetDecayRCut(5.);
    taskHighPtDeDx->SetMinPtV0(1.0);
    taskHighPtDeDx->SetMassCut(0.1); 
    taskHighPtDeDx->SetRejectKinks(kTRUE);
    taskHighPtDeDx->SetSigmaDedxCut(kFALSE);
    
    //Set Filters
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
        
    mgr->AddTask(taskHighPtDeDx);
    
  }
  
if(typerun==3){//pp analysis
  
  AliAnalysisTaskHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx("taskHighPtDeDx");

    //Set analysis details
    taskHighPtDeDx->SetAnalysisType("AOD");
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kFALSE); //for pbpb and ppb
    taskHighPtDeDx->SetAnalysisRun2(kTRUE);
    taskHighPtDeDx->SetTrigger1(kTriggerInt); 
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetDebugLevel(0);
    
    //Set event details
    taskHighPtDeDx->SetCentFrameworkAliCen(CentFrameworkAliCen); //kTRUE:AliCentrality, kFALSE:AliMultSelection
    taskHighPtDeDx->SetCentDetector("V0M");
    taskHighPtDeDx->SetMinCent(minc);
    taskHighPtDeDx->SetMaxCent(maxc);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetContributorsVtxCut(0);
    taskHighPtDeDx->SetContributorsVtxSPDCut(0);
    taskHighPtDeDx->SetPileupCut(999.); //Correlation between global Zvtx and SPD Zvtx: use 999. for no cut
    taskHighPtDeDx->SetVtxR2Cut(10.); //use 10. for no cut
    
    //Set trigger particle tracks, v0s and daughter track details
    taskHighPtDeDx->SetMinPt(5.0); //trigger tracks
    taskHighPtDeDx->SetLowPtFraction(0.0); //keep x.x% of tracks below min pt
    taskHighPtDeDx->SetEtaCut(0.8); 
    taskHighPtDeDx->SetCrossedRowsCut(70.); //use 0. for no cut
    taskHighPtDeDx->SetCrossedOverFindableCut(0.8); //use 0. for no cut
    taskHighPtDeDx->SetCosPACut(0.97);
    taskHighPtDeDx->SetDecayRCut(0.5);
    taskHighPtDeDx->SetMinPtV0(1.0);
    taskHighPtDeDx->SetMassCut(0.1); 
    taskHighPtDeDx->SetRejectKinks(kTRUE);
    taskHighPtDeDx->SetSigmaDedxCut(kFALSE);
    
    //Set Filters
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);

    mgr->AddTask(taskHighPtDeDx);
  
  
 }


if(typerun==4){//ppb analysis
  cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<    Runing pPb    <<<<<<<<<<<<<<"<<endl;
    
  AliAnalysisTaskHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx("taskHighPtDeDx");

    //Set analysis details
    taskHighPtDeDx->SetAnalysisType("AOD");
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE); //for pbpb and ppb
    taskHighPtDeDx->SetAnalysisRun2(kTRUE);
    taskHighPtDeDx->SetTrigger1(kTriggerInt); 
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetDebugLevel(0);
    
    //Set event details
    taskHighPtDeDx->SetCentFrameworkAliCen(CentFrameworkAliCen); //kTRUE:AliCentrality, kFALSE:AliMultSelection
    taskHighPtDeDx->SetCentDetector("V0A");
    taskHighPtDeDx->SetMinCent(minc);
    taskHighPtDeDx->SetMaxCent(maxc);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetContributorsVtxCut(0);
    taskHighPtDeDx->SetContributorsVtxSPDCut(0);
    taskHighPtDeDx->SetPileupCut(999.); //Correlation between global Zvtx and SPD Zvtx: use 999. for no cut
    taskHighPtDeDx->SetVtxR2Cut(10.); //use 10. for no cut
    
    //Set trigger particle tracks, v0s and daughter track details
    taskHighPtDeDx->SetMinPt(5.0); //trigger tracks
    taskHighPtDeDx->SetLowPtFraction(0.0); //keep x.x% of tracks below min pt
    taskHighPtDeDx->SetEtaCut(0.8); 
    taskHighPtDeDx->SetCrossedRowsCut(70.); //use 0. for no cut
    taskHighPtDeDx->SetCrossedOverFindableCut(0.8); //use 0. for no cut
    taskHighPtDeDx->SetCosPACut(0.97);
    taskHighPtDeDx->SetDecayRCut(0.5);
    taskHighPtDeDx->SetMinPtV0(1.0);
    taskHighPtDeDx->SetMassCut(0.1); 
    taskHighPtDeDx->SetRejectKinks(kTRUE);
    taskHighPtDeDx->SetSigmaDedxCut(kFALSE);

    //Set Filters
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    
    mgr->AddTask(taskHighPtDeDx);
    
 }


// Create ONLY the output containers for the data produced by the
// task.  Get and connect other common input/output containers via
// the manager as below
//=======================================================================
 
 TString outputFileName = Form("%s", AliAnalysisManager::GetCommonFileName());
 AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
 
 mgr->ConnectInput(taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
 mgr->ConnectOutput(taskHighPtDeDx, 1, cout_hist);
 
 
 // Return task pointer at the end
 return taskHighPtDeDx;
 
}
