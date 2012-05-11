AliAnalysisTask* AddTask(Bool_t AnalysisMC, const Char_t* taskname, Int_t typerun, UInt_t kTriggerInt[], Float_t minc[], Float_t maxc[] )
{
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access method
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
  AliESDtrackCuts* esdTrackCutsGolden = 
    AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
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
    
    
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx[6];
    for( Int_t i=0; i<6; ++i ){
      taskHighPtDeDx[i]=0;
      Char_t TaskName[256]={0};
      sprintf(TaskName,"%s_%1.0f_%1.0f",taskname,minc[i],maxc[i]);
      
      taskHighPtDeDx[i] = new AliAnalysisTaskHighPtDeDx(TaskName);
      TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
      taskHighPtDeDx[i]->SetAnalysisType(type);
      taskHighPtDeDx[i]->SetAnalysisMC(AnalysisMC);
      taskHighPtDeDx[i]->SetAnalysisPbPb(kTRUE);
      taskHighPtDeDx[i]->SetProduceTPCBranch(kFALSE);
      taskHighPtDeDx[i]->SetDebugLevel(0);
      taskHighPtDeDx[i]->SetEtaCut(0.8);
      taskHighPtDeDx[i]->SetVtxCut(10.0);
      taskHighPtDeDx[i]->SetMinPt(0.0); // default 2.0
      taskHighPtDeDx[i]->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
      taskHighPtDeDx[i]->SetTreeOption(1);
      taskHighPtDeDx[i]->SetTrigger1(kTriggerInt[0]);
      taskHighPtDeDx[i]->SetTrigger2(kTriggerInt[1]);
      taskHighPtDeDx[i]->SetMinCent(minc[i]);
      taskHighPtDeDx[i]->SetMaxCent(maxc[i]);
      //Set Filtesr
      taskHighPtDeDx[i]->SetTrackFilterGolden(trackFilterGolden);
      taskHighPtDeDx[i]->SetTrackFilter(trackFilter0);
      taskHighPtDeDx[i]->SetTrackFilterTPC(trackFilterTPC);
      
      mgr->AddTask(taskHighPtDeDx[i]);
      
    }
    
    // Create ONLY the output containers for the data produced by the
    // task.  Get and connect other common input/output containers via
    // the manager as below
    //=======================================================================
    AliAnalysisDataContainer *cout_hist[6];
    for( Int_t i=0; i<6; ++i ){
      
      cout_hist[i]=0;
      Char_t outFileName[256]={0};
      sprintf(outFileName,"%s_Tree_%1.0f_%1.0f.root",taskname,minc[i],maxc[i]);
      //AliAnalysisDataContainer *cout_hist    = 0;
      
      cout_hist[i] = mgr->CreateContainer(Form("output_%1.0f_%1.0f",minc[i],maxc[i]), TList::Class(), AliAnalysisManager::kOutputContainer, outFileName);
      mgr->ConnectInput (taskHighPtDeDx[i], 0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskHighPtDeDx[i], 1, cout_hist[i]);
      
    }  
    
    // Return task pointer at the end
    return taskHighPtDeDx[0];
  }
  if(typerun==3){//pp analysis
  
    
    AliAnalysisTaskHighPtDeDx* taskHighPtDeDx = new AliAnalysisTaskHighPtDeDx("taskHighPtDeDxpp");
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    taskHighPtDeDx->SetAnalysisType(type);
    taskHighPtDeDx->SetAnalysisMC(AnalysisMC);
    taskHighPtDeDx->SetAnalysisPbPb(kFALSE);
    taskHighPtDeDx->SetProduceTPCBranch(kFALSE);
    taskHighPtDeDx->SetDebugLevel(0);
    taskHighPtDeDx->SetEtaCut(0.8);
    taskHighPtDeDx->SetVtxCut(10.0);
    taskHighPtDeDx->SetMinPt(0.0); // default 2.0
    taskHighPtDeDx->SetLowPtFraction(0.01); // keep 1% of tracks below min pt
    taskHighPtDeDx->SetTreeOption(1);
    taskHighPtDeDx->SetTrigger1(kTriggerInt[0]);
    taskHighPtDeDx->SetTrigger2(kTriggerInt[1]);
    //Set Filtesr
    taskHighPtDeDx->SetTrackFilterGolden(trackFilterGolden);
    taskHighPtDeDx->SetTrackFilter(trackFilter0);
    taskHighPtDeDx->SetTrackFilterTPC(trackFilterTPC);
    
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
  
  
}
