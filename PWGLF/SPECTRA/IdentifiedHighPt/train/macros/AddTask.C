/*
  Last update (new version): 
  20 aug 2015: clean up
  22 aug 2015: match pre defined output

*/

AliAnalysisTask* AddTask(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun, UInt_t kTriggerInt, Float_t minc, Float_t maxc)
{
   
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access methodh
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
  /*
  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = 
    AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);
  */

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
  Bool_t clusterCut=0;
  Bool_t selPrimaries=kTRUE;
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
  esdTrackCutsL0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  if(selPrimaries) {
    // 7*(0.0026+0.0050/pt^1.01)
    esdTrackCutsL0->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  }
  esdTrackCutsL0->SetMaxDCAToVertexZ(2);
  esdTrackCutsL0->SetDCAToVertex2D(kFALSE);
  esdTrackCutsL0->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsL0->SetMaxChi2PerClusterITS(1e10);
  esdTrackCutsL0->SetMaxChi2TPCConstrainedGlobal(1e10);


  trackFilter0->AddCuts(esdTrackCutsL0);

  
  
  AliAnalysisFilter* trackFilterTPC = new AliAnalysisFilter("trackFilterTPC");
  AliESDtrackCuts* esdTrackCutsTPC = 
    AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackFilterTPC->AddCuts(esdTrackCutsTPC);



  // Create the task and configure it 
  //========================================================================
  if(typerun==2){//heavy ion analysis
    
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx;
    
    taskHighPtDeDx=0;
    Char_t TaskName[256]={0};
    sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc,maxc);
    
    taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx(TaskName);
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetMinPt(0.0); 
    taskHighPtDeDx->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
    taskHighPtDeDx->SetTreeOption(1);
    taskHighPtDeDx->SetTrigger1(kTriggerInt); 
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    // taskHighPtDeDx->SetTrigger1(AliVEvent::kMB); 
    // taskHighPtDeDx->SetTrigger2(AliVEvent::kMB);
    taskHighPtDeDx->SetMinCent(minc);
    taskHighPtDeDx->SetMaxCent(maxc);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    
    //V0's cuts
    taskHighPtDeDx->SetMassCut(0.1);     
    taskHighPtDeDx->SetRequireRecV0(kTRUE); // def: kTRUE
    //taskHighPtDeDx->SetStoreMcIn(kFALSE);     // def: kFALSE
    taskHighPtDeDx->SetStoreMcIn(kTRUE);     // def: kFALSE
    
    mgr->AddTask(taskHighPtDeDx);
    
  }
  
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //=======================================================================
  
    TString outputFileName = Form("%s", AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
 
  mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHighPtDeDx, 1, cout_hist);


 // Return task pointer at the end
  return taskHighPtDeDx;
}
if(typerun==3){//pp analysis
  
  
  AliAnalysisTaskHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx("taskHighPtDeDxpp");
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  taskHighPtDeDx->SetAnalysisType(type);
  taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kFALSE);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetMinPt(0.0); // default 2.0
    taskHighPtDeDx->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
    taskHighPtDeDx->SetTreeOption(1);
    taskHighPtDeDx->SetTrigger1(kTriggerInt);
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    // taskHighPtDeDx->SetTrigger1(AliVEvent::kMB);
    // taskHighPtDeDx->SetTrigger2(AliVEvent::kMB);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    //V0's cuts
    taskHighPtDeDx->SetMassCut(0.1);         // def: 0.03
    taskHighPtDeDx->SetRequireRecV0(kTRUE); // def: kTRUE
    taskHighPtDeDx->SetStoreMcIn(kFALSE);     // def: kFALSE

    mgr->AddTask(taskHighPtDeDx);
      
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
  
    AliAnalysisDataContainer *cout_histdedx;

    cout_histdedx=0;
    Char_t outFileName[256]={0};
    sprintf(outFileName,"%s_Tree.root",taskname);
    //AliAnalysisDataContainer *cout_hist    = 0;
 
    cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDx, 1, cout_histdedx);
 
    // Return task pointer at the end
    return taskHighPtDeDx;
    

    
    
    
  }



 if(typerun==4){//pPb analysis
   cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<                 Runing pPb    <<<<<<<<<<<<<<"<<endl;
    
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx("taskHighPtDeDxpp");
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kTRUE);
    taskHighPtDeDx->SetMinCent(-200);
    taskHighPtDeDx->SetMaxCent(200);
    taskHighPtDeDx->SetProduceVZEROBranch(kTRUE);
    taskHighPtDeDx->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetMinPt(0.0); // default 2.0
    taskHighPtDeDx->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
    taskHighPtDeDx->SetTreeOption(1);
    taskHighPtDeDx->SetTrigger1(kTriggerInt);
    taskHighPtDeDx->SetTrigger2(kTriggerInt);
    // taskHighPtDeDx->SetTrigger1(AliVEvent::kMB);
    // taskHighPtDeDx->SetTrigger2(AliVEvent::kMB);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    //V0's cuts
    taskHighPtDeDx->SetMassCut(0.1);       
    taskHighPtDeDx->SetRequireRecV0(kTRUE); // def: kTRUE
    taskHighPtDeDx->SetStoreMcIn(kTRUE);     // def: kFALSE

    mgr->AddTask(taskHighPtDeDx);
      
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
  
    AliAnalysisDataContainer *cout_histdedx;

    cout_histdedx=0;
    Char_t outFileName[256]={0};
    sprintf(outFileName,"%s_Tree.root",taskname);
    //AliAnalysisDataContainer *cout_hist    = 0;
 
    cout_histdedx = mgr->CreateContainer("outputdedx", TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
    mgr->ConnectInput (taskHighPtDeDx, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHighPtDeDx, 1, cout_histdedx);
 
    // Return task pointer at the end
    return taskHighPtDeDx;
    

    
    
    
 }
